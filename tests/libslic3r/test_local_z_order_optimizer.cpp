#include <catch2/catch.hpp>
#include "libslic3r/LocalZOrderOptimizer.hpp"

using namespace Slic3r::LocalZOrderOptimizer;

// ---------------------------------------------------------------------------
// bucket_contains_extruder
// ---------------------------------------------------------------------------
TEST_CASE("bucket_contains_extruder basics", "[LocalZOrderOptimizer]")
{
    std::vector<unsigned int> bucket = {1, 3, 5};
    CHECK(bucket_contains_extruder(bucket, 1));
    CHECK(bucket_contains_extruder(bucket, 3));
    CHECK(bucket_contains_extruder(bucket, 5));
    CHECK_FALSE(bucket_contains_extruder(bucket, 2));
    CHECK_FALSE(bucket_contains_extruder(bucket, 0));
    CHECK_FALSE(bucket_contains_extruder(bucket, -1));
    CHECK_FALSE(bucket_contains_extruder({}, 1));
}

// ---------------------------------------------------------------------------
// order_bucket_extruders
// ---------------------------------------------------------------------------
TEST_CASE("order_bucket_extruders rotates current to front",
          "[LocalZOrderOptimizer]")
{
    std::vector<unsigned int> ext = {1, 2, 3};
    auto result = order_bucket_extruders(ext, 2);
    REQUIRE(result.size() == 3);
    CHECK(result.front() == 2);
}

TEST_CASE("order_bucket_extruders puts preferred_last at end",
          "[LocalZOrderOptimizer]")
{
    std::vector<unsigned int> ext = {1, 2, 3, 4};
    auto result = order_bucket_extruders(ext, 1, 3);
    REQUIRE(result.size() == 4);
    CHECK(result.front() == 1);
    CHECK(result.back() == 3);
}

TEST_CASE("order_bucket_extruders handles empty input",
          "[LocalZOrderOptimizer]")
{
    auto result = order_bucket_extruders({}, 1);
    CHECK(result.empty());
}

TEST_CASE("order_bucket_extruders deduplicates", "[LocalZOrderOptimizer]")
{
    std::vector<unsigned int> ext = {1, 1, 2, 2, 3};
    auto result = order_bucket_extruders(ext, 1);
    CHECK(result.size() <= 3); // unique removes adjacent duplicates
}

// ---------------------------------------------------------------------------
// order_pass_group
// ---------------------------------------------------------------------------
TEST_CASE("order_pass_group prioritizes groups with current extruder",
          "[LocalZOrderOptimizer]")
{
    std::vector<std::vector<unsigned int>> groups = {
        {2, 3},
        {1, 2},
        {3, 4}
    };

    auto result = order_pass_group(groups, 1);
    REQUIRE(result.size() == 3);
    // Group containing extruder 1 should come first
    CHECK(result[0] == 1); // groups[1] = {1,2}
}

TEST_CASE("order_pass_group produces valid permutation",
          "[LocalZOrderOptimizer]")
{
    std::vector<std::vector<unsigned int>> groups = {
        {1},
        {2},
        {3},
        {4}
    };

    auto result = order_pass_group(groups, 2);
    REQUIRE(result.size() == 4);

    // Should be a permutation of {0,1,2,3}
    std::vector<size_t> sorted = result;
    std::sort(sorted.begin(), sorted.end());
    CHECK(sorted == std::vector<size_t>{0, 1, 2, 3});

    // Group with extruder 2 first
    CHECK(result[0] == 1);
}

TEST_CASE("order_pass_group with single group", "[LocalZOrderOptimizer]")
{
    std::vector<std::vector<unsigned int>> groups = {{1, 2, 3}};
    auto result = order_pass_group(groups, 1);
    REQUIRE(result.size() == 1);
    CHECK(result[0] == 0);
}

TEST_CASE("order_pass_group handles empty groups", "[LocalZOrderOptimizer]")
{
    std::vector<std::vector<unsigned int>> groups;
    auto result = order_pass_group(groups, 1);
    CHECK(result.empty());
}

// ---------------------------------------------------------------------------
// Deterministic ordering: same input always produces same output
// ---------------------------------------------------------------------------
TEST_CASE("order_pass_group is deterministic", "[LocalZOrderOptimizer]")
{
    std::vector<std::vector<unsigned int>> groups = {
        {1, 3},
        {2, 4},
        {1, 4},
        {3, 2}
    };

    auto r1 = order_pass_group(groups, 1);
    auto r2 = order_pass_group(groups, 1);
    CHECK(r1 == r2);
}
