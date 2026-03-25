#include <catch2/catch.hpp>
#include "libslic3r/MixedFilament.hpp"

using namespace Slic3r;

// ---------------------------------------------------------------------------
// MixedFilamentManager – auto_generate
// ---------------------------------------------------------------------------
TEST_CASE("MixedFilamentManager auto_generate creates pairwise combinations",
          "[MixedFilament]")
{
    MixedFilamentManager mgr;
    std::vector<std::string> colors = {"#FF0000", "#00FF00", "#0000FF"};
    mgr.auto_generate(colors);

    const auto &mixed = mgr.mixed_filaments();
    REQUIRE(mixed.size() == 3); // C(3,2) = 3

    // Each pair should reference valid 1-based IDs
    REQUIRE(mixed[0].component_a == 1);
    REQUIRE(mixed[0].component_b == 2);
    REQUIRE(mixed[1].component_a == 1);
    REQUIRE(mixed[1].component_b == 3);
    REQUIRE(mixed[2].component_a == 2);
    REQUIRE(mixed[2].component_b == 3);

    // All enabled by default
    for (const auto &mf : mixed) {
        CHECK(mf.enabled);
        CHECK_FALSE(mf.deleted);
        CHECK_FALSE(mf.custom);
        CHECK(mf.origin_auto);
        CHECK(mf.stable_id != 0);
    }
}

TEST_CASE("auto_generate with 2 filaments produces 1 pair", "[MixedFilament]")
{
    MixedFilamentManager mgr;
    mgr.auto_generate({"#FF0000", "#00FF00"});
    REQUIRE(mgr.mixed_filaments().size() == 1);
    REQUIRE(mgr.enabled_count() == 1);
    REQUIRE(mgr.total_filaments(2) == 3);
}

TEST_CASE("auto_generate with 1 filament produces nothing", "[MixedFilament]")
{
    MixedFilamentManager mgr;
    mgr.auto_generate({"#FF0000"});
    REQUIRE(mgr.mixed_filaments().empty());
    REQUIRE(mgr.total_filaments(1) == 1);
}

TEST_CASE("auto_generate with 4 filaments produces C(4,2)=6 pairs",
          "[MixedFilament]")
{
    MixedFilamentManager mgr;
    mgr.auto_generate({"#FF0000", "#00FF00", "#0000FF", "#FFFFFF"});
    REQUIRE(mgr.mixed_filaments().size() == 6);
}

// ---------------------------------------------------------------------------
// Nozzle-diameter-aware pairing
// ---------------------------------------------------------------------------
TEST_CASE("auto_generate filters by nozzle diameter", "[MixedFilament]")
{
    MixedFilamentManager mgr;
    // 3 filaments on 0.4mm, 1 on 0.6mm
    std::vector<std::string> colors = {"#FF0000", "#00FF00", "#0000FF", "#FFFFFF"};
    std::vector<double> diameters   = {0.4, 0.4, 0.4, 0.6};
    mgr.auto_generate(colors, diameters);

    // Only C(3,2)=3 pairs among the 0.4mm group, 0 from the lone 0.6mm
    REQUIRE(mgr.mixed_filaments().size() == 3);
    for (const auto &mf : mgr.mixed_filaments()) {
        CHECK(mf.component_a <= 3);
        CHECK(mf.component_b <= 3);
    }
}

TEST_CASE("auto_generate all same diameter pairs all", "[MixedFilament]")
{
    MixedFilamentManager mgr;
    std::vector<std::string> colors = {"#FF0000", "#00FF00", "#0000FF"};
    std::vector<double> diameters   = {0.4, 0.4, 0.4};
    mgr.auto_generate(colors, diameters);
    REQUIRE(mgr.mixed_filaments().size() == 3); // C(3,2) = 3
}

TEST_CASE("auto_generate with two diameter groups", "[MixedFilament]")
{
    MixedFilamentManager mgr;
    // 2 filaments on 0.4 + 2 filaments on 0.6
    std::vector<std::string> colors = {"#A", "#B", "#C", "#D"};
    std::vector<double> diameters   = {0.4, 0.4, 0.6, 0.6};
    mgr.auto_generate(colors, diameters);

    // C(2,2)=1 for 0.4 group + C(2,2)=1 for 0.6 group = 2
    REQUIRE(mgr.mixed_filaments().size() == 2);
    CHECK(mgr.mixed_filaments()[0].component_a == 1);
    CHECK(mgr.mixed_filaments()[0].component_b == 2);
    CHECK(mgr.mixed_filaments()[1].component_a == 3);
    CHECK(mgr.mixed_filaments()[1].component_b == 4);
}

// ---------------------------------------------------------------------------
// Stable ID preservation across regenerations
// ---------------------------------------------------------------------------
TEST_CASE("auto_generate preserves stable_id and enabled state",
          "[MixedFilament]")
{
    MixedFilamentManager mgr;
    mgr.auto_generate({"#FF0000", "#00FF00", "#0000FF"});

    // Disable pair (1,3)
    auto &mixed = mgr.mixed_filaments();
    REQUIRE(mixed.size() == 3);
    uint64_t id_pair_1_3 = mixed[1].stable_id;
    mixed[1].enabled = false;

    // Re-generate (e.g. after filament colour change)
    mgr.auto_generate({"#FF0000", "#00FF00", "#0000FF"});

    const auto &after = mgr.mixed_filaments();
    REQUIRE(after.size() == 3);
    // Pair (1,3) should still be disabled and retain its stable_id
    CHECK(after[1].component_a == 1);
    CHECK(after[1].component_b == 3);
    CHECK(after[1].stable_id == id_pair_1_3);
    CHECK_FALSE(after[1].enabled);
}

// ---------------------------------------------------------------------------
// Resolve
// ---------------------------------------------------------------------------
TEST_CASE("Resolve returns physical ID for non-mixed filaments",
          "[MixedFilament]")
{
    MixedFilamentManager mgr;
    mgr.auto_generate({"#FF0000", "#00FF00"});

    // Filament 1 (physical) should pass through
    CHECK(mgr.resolve(1, 2, 0) == 1);
    CHECK(mgr.resolve(2, 2, 0) == 2);
}

TEST_CASE("Resolve alternates components for mixed filament",
          "[MixedFilament]")
{
    MixedFilamentManager mgr;
    mgr.auto_generate({"#FF0000", "#00FF00"});

    // Mixed filament is ID 3 (num_physical=2, first virtual)
    unsigned int layer0 = mgr.resolve(3, 2, 0);
    unsigned int layer1 = mgr.resolve(3, 2, 1);
    // With ratio_a=1, ratio_b=1, should alternate between 1 and 2
    CHECK(layer0 != layer1);
    CHECK((layer0 == 1 || layer0 == 2));
    CHECK((layer1 == 1 || layer1 == 2));
}

// ---------------------------------------------------------------------------
// is_mixed
// ---------------------------------------------------------------------------
TEST_CASE("is_mixed correctly identifies virtual filament IDs",
          "[MixedFilament]")
{
    MixedFilamentManager mgr;
    mgr.auto_generate({"#FF0000", "#00FF00", "#0000FF"});

    // 3 physical, 3 mixed → IDs 4, 5, 6 are mixed
    CHECK_FALSE(mgr.is_mixed(1, 3));
    CHECK_FALSE(mgr.is_mixed(2, 3));
    CHECK_FALSE(mgr.is_mixed(3, 3));
    CHECK(mgr.is_mixed(4, 3));
    CHECK(mgr.is_mixed(5, 3));
    CHECK(mgr.is_mixed(6, 3));
    CHECK_FALSE(mgr.is_mixed(7, 3)); // out of range
}

// ---------------------------------------------------------------------------
// total_filaments and enabled_count
// ---------------------------------------------------------------------------
TEST_CASE("total_filaments reflects enabled count", "[MixedFilament]")
{
    MixedFilamentManager mgr;
    mgr.auto_generate({"#FF0000", "#00FF00", "#0000FF"});

    REQUIRE(mgr.enabled_count() == 3);
    REQUIRE(mgr.total_filaments(3) == 6);

    // Disable one
    mgr.mixed_filaments()[0].enabled = false;
    CHECK(mgr.enabled_count() == 2);
    CHECK(mgr.total_filaments(3) == 5);
}

// ---------------------------------------------------------------------------
// remove_physical_filament
// ---------------------------------------------------------------------------
TEST_CASE("remove_physical_filament removes affected pairs and shifts IDs",
          "[MixedFilament]")
{
    MixedFilamentManager mgr;
    mgr.auto_generate({"#FF0000", "#00FF00", "#0000FF"});
    REQUIRE(mgr.mixed_filaments().size() == 3);

    // Remove filament 2 → pairs containing 2 removed, IDs shifted
    mgr.remove_physical_filament(2);

    const auto &after = mgr.mixed_filaments();
    // Only pair (1,3) survives, renumbered to (1,2) after shift
    REQUIRE(after.size() == 1);
    CHECK(after[0].component_a == 1);
    CHECK(after[0].component_b == 2);
}

// ---------------------------------------------------------------------------
// Custom entries and serialization
// ---------------------------------------------------------------------------
TEST_CASE("Custom entry round-trips through serialize/load",
          "[MixedFilament]")
{
    MixedFilamentManager mgr;
    std::vector<std::string> colors = {"#FF0000", "#00FF00", "#0000FF"};
    mgr.auto_generate(colors);
    mgr.add_custom_filament(1, 3, 30, colors);

    std::string serialized = mgr.serialize_custom_entries();
    CHECK_FALSE(serialized.empty());

    MixedFilamentManager mgr2;
    mgr2.auto_generate(colors);
    mgr2.load_custom_entries(serialized, colors);

    // Should have the same total count
    CHECK(mgr2.mixed_filaments().size() == mgr.mixed_filaments().size());
}

// ---------------------------------------------------------------------------
// blend_color
// ---------------------------------------------------------------------------
TEST_CASE("blend_color produces valid hex output", "[MixedFilament]")
{
    std::string result = MixedFilamentManager::blend_color("#FF0000", "#0000FF", 1, 1);
    // Should be a valid #RRGGBB string
    REQUIRE(result.size() == 7);
    CHECK(result[0] == '#');
}

// ---------------------------------------------------------------------------
// normalize_manual_pattern
// ---------------------------------------------------------------------------
TEST_CASE("normalize_manual_pattern cleans input", "[MixedFilament]")
{
    CHECK(MixedFilamentManager::normalize_manual_pattern("1,2,1,2") == "1212");
    CHECK(MixedFilamentManager::normalize_manual_pattern("AABB") == "1122");
    CHECK(MixedFilamentManager::normalize_manual_pattern("") == "");
}

TEST_CASE("mix_percent_from_manual_pattern computes B-ratio",
          "[MixedFilament]")
{
    // "1122" = 50% B
    CHECK(MixedFilamentManager::mix_percent_from_manual_pattern("1122") == 50);
    // "111222" = 50% B .. though this depends on impl
    CHECK(MixedFilamentManager::mix_percent_from_manual_pattern("12") == 50);
}
