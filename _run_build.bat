@echo off
REM Remove Strawberry's c\bin from PATH (has cmake 3.29 that conflicts with VS cmake)
set "PATH=%PATH:C:\Strawberry\c\bin;=%"
set "PATH=%PATH:;C:\Strawberry\c\bin=%"
set "PATH=%PATH:C:\Strawberry\perl\site\bin;=%"
set "PATH=%PATH:;C:\Strawberry\perl\site\bin=%"
set "PATH=%PATH:C:\Strawberry\perl\bin;=%"
set "PATH=%PATH:;C:\Strawberry\perl\bin=%"
REM Add only Strawberry perl (not cmake) back to PATH
set "PATH=C:\Strawberry\perl\bin;%PATH%"
REM Add deps bin to PATH for pkg-config
set "DEPS_DIR=a:\Users\Iwo\Nextcloud\git\H2C-FullSpectrum\BambuStudio\deps\build\out_deps\usr\local"
set "PATH=%DEPS_DIR%\bin;%PATH%"
REM Set PKG_CONFIG_PATH so pkg-config finds the .pc files from deps
set "PKG_CONFIG_PATH=%DEPS_DIR%\lib\pkgconfig"
cd /d "a:\Users\Iwo\Nextcloud\git\H2C-FullSpectrum\BambuStudio\build"
REM Cache already clean from previous run, skip deletion
cd /d "a:\Users\Iwo\Nextcloud\git\H2C-FullSpectrum\BambuStudio"
call build_win.bat -s app-dirty -d deps\build\out_deps -c Release
