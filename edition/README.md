# Vector-based terrain modelling - Edition

## Requirements

This project has been tested on Windows with CLion 2025 and depends on C++20, CMake, Qt 6.9.1. While it may work on other platform, it has not been tested for it.

## Installation

We rely on CMake for its configuration. Hence, executable can be built on Windows with:
```
git clone --recursive https://github.com/simonperche/vector-based_terrain_modelling.git
cd vector-based_terrain_modelling/edition
mkdir build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config "Release"
```

Since this project uses Qt, you must ensure that the required Qt dynamic libraries are available in the output directory.
On Windows, use `windeployqt` (bundled with Qt) to automatically copy the necessary DLLs:
```
{QT_PATH}\{version}\{compiler}\bin\windeployqt.exe Release\VectorTerrains.exe
```

## Usage

_Coming soon..._