# Vector-based terrain modelling - Edition

## Requirements

This project has been tested on Windows with CLion 2025 and depends on C++20, CMake, Qt 6.9.1. While it may work on other platform, it has not been tested for it.

### On linux
You need to install the following packages
```qt6-base-dev qt6-declarative-dev qt6-tools-dev qt6-charts-dev qt6-opengl-dev libglew-dev libxkbcommon-dev build-essential cmake```

## Installation

We rely on CMake for its configuration. Hence, executable can be built on Windows with:

### Windows
```
git clone --recursive https://github.com/simonperche/vector-based_terrain_modelling.git
cd vector-based_terrain_modelling/edition
mkdir build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config "Release"
```

### Linux
```
git clone --recursive https://github.com/simonperche/vector-based_terrain_modelling.git
cd vector-based_terrain_modelling/edition
mkdir build
cmake -S . -B build
cmake --build build -j$(nproc)
```

Since this project uses Qt, you must ensure that the required Qt dynamic libraries are available in the output directory.
On Windows, use `windeployqt` (bundled with Qt) to automatically copy the necessary DLLs:
```
{QT_PATH}\{version}\{compiler}\bin\windeployqt.exe Release\VectorTerrains.exe
```

Linux will handle it automatically.


## Usage

### Global shortcuts

- Middle click : rotate the camera (Use Ctrl+MMB when using Crest Graph)
- Shift + middle click : pan the camera
- Left click : use the current tool

#### Tool shortcuts

##### Graph tool

- Mouse click on a vertex : move the selected subgraph
- Alt + wheel : extend the depth in the graph
- Ctrl + wheel : change the amplitude of the selected subgraph
- Ctrl + alt + wheel : change the influence region (if enabled)

##### Local modifications

- Ctrl + wheel : scale the area
- Amplitude : click on the terrain + move up and down

##### Region

- Click and drag : create a region
- Ctrl + wheel : rotate the region
- Ctrl + alt + wheel : scale the region

### UI

![](../docs/edition_ui.jpg)

1. Input/Output : `Load primitives` loads [.npy | .csv] files. `Save primitives` saves a .csv file. `Add primitives` loads a ground truth heightmap and add details primitives to the terrain. `Export high-res` exports the current terrain as a heightmap, with the possibility to load a ground truth heightmaps for the details primitives.
2. View parameters : `Nb primitivies` shows X primitives. `Show influence of the primitives` displays an ellipse around each primitive at 3$\sigma$. `Amplitude` and `Details level` change the amplitude of respectively the entire terrain or only the details primitives.
3. Tools : manipulate the terrain (parameters of each tool will display in 3.2 when selected)
   - `Local modifications` : erase/change amplitude in a local area.
   - `Region` : select and move a portion of terrain, which can be saved as a `.csv`, and loaded as a template.
   - `Crest graph` : geomorphology-based graph to move and edit structures
4. Miscellaneous : `Save logs` saves each modification into `./logs` folder, checked by default. `Load Heightfield` loads a .png file for display into the app. `Reload shader` reloads each `.glsl` files from the disk, useful for debugging.
