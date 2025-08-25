# in order to propoer work, you have to set this env variable to point to your Qt install dir
# CMAKE_PREFIX_PATH
cmake_minimum_required(VERSION 3.16)
set(TARGET_NAME VectorTerrains)
project(${TARGET_NAME} VERSION 1.0.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_FLAGS "-mavx -mavx2 -march=native")
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(Qt6 REQUIRED COMPONENTS Core)
find_package(Qt6 REQUIRED COMPONENTS Gui)
find_package(Qt6 REQUIRED COMPONENTS Widgets)
find_package(Qt6 REQUIRED COMPONENTS Charts)
find_package(Qt6 REQUIRED COMPONENTS PrintSupport)

#set(CMAKE_AUTOUIC ON)

# OpenGL & GLEW library
find_package(OpenGL)
find_package(GLEW)
if (NOT GLEW_FOUND)
    message("Adding GLEW manually...")
    set(PLATFORM x64)
    set(GLEW_INCLUDE_DIRS $ENV{ARCHESLIBDIR}/LibPlatform/LibGL/$(PLATFORM)/Include)
    set(GLEW_LIBRARIES $ENV{ARCHESLIBDIR}/LibPlatform/LibGL/$(PLATFORM)/Library/glew32.lib)
    set(GLEW_FOUND true)
endif ()

set(HEADER_FILES
        Include/GaussianTerrainRaytracingWidget.h
        Include/graph.h
        Include/Kernels.h
        Include/MayaSimpleRendererColors.h
        Include/npy.hpp
        Include/qte.h
        Include/utils.h

        Include/Tools/Tool.h
        Include/Tools/ToolBrush.h
        Include/Tools/ToolCurve.h
        Include/Tools/ToolEdit.h
        Include/Tools/ToolGraph.h

        Include/Kernels/FactoryKernel.h
        Include/Kernels/DetailsKernel.h
        Include/Kernels/GaussianKernel.h
        Include/Kernels/Kernel.h
)

set(SOURCE_FILES
        Source/GaussianTerrainRaytracingWidget.cpp
        Source/graph.cpp
        Source/Kernels.cpp
        Source/main.cpp
        Source/MayaSimpleRendererColors.cpp
        Source/qte-main-window.cpp
        Source/utils.cpp

        Source/Tools/ToolBrush.cpp
        Source/Tools/ToolCurve.cpp
        Source/Tools/ToolEdit.cpp
        Source/Tools/ToolGraph.cpp

        Source/Kernels/DetailsKernel.cpp
        Source/Kernels/GaussianKernel.cpp
        Source/Kernels/Kernel.cpp
)

qt_wrap_ui(SOURCE_FILES UI/main.ui)
qt_wrap_cpp(SOURCE_FILES Include/qte.h)
qt_wrap_cpp(SOURCE_FILES Include/GaussianTerrainRaytracingWidget.h)

qt_add_executable(${TARGET_NAME} ${SOURCE_FILES})

target_compile_definitions(${TARGET_NAME} PUBLIC SOLUTION_DIR="${CMAKE_CURRENT_SOURCE_DIR}/")
target_include_directories(${TARGET_NAME} PRIVATE Include)
target_include_directories(${TARGET_NAME} PRIVATE $ENV{ARCHESLIBDIR}/LibCore/LibCore/Include)
target_include_directories(${TARGET_NAME} PRIVATE $ENV{ARCHESLIBDIR}/LibMaya/LibMaya/Include)
target_include_directories(${TARGET_NAME} PRIVATE $ENV{ARCHESLIBDIR}/LibAssets/LibAssets/Include)
target_include_directories(${TARGET_NAME} PRIVATE $ENV{ARCHESLIBDIR}/LibHeightfield/LibHeightfield/Include)

target_include_directories(${TARGET_NAME} PRIVATE ${GLEW_INCLUDE_DIRS})

# this assumes that the build dir is "build" and the libraries have been built
target_link_directories(${TARGET_NAME} BEFORE PUBLIC $ENV{ARCHESLIBDIR}/LibAssets/build)
target_link_directories(${TARGET_NAME} BEFORE PUBLIC $ENV{ARCHESLIBDIR}/LibMaya/build)
target_link_directories(${TARGET_NAME} BEFORE PUBLIC $ENV{ARCHESLIBDIR}/LibHeightfield/build)
target_link_directories(${TARGET_NAME} BEFORE PUBLIC $ENV{ARCHESLIBDIR}/LibCore/build)
target_link_directories(${TARGET_NAME} BEFORE PUBLIC $ENV{GLEW_INCLUDE_DIRS})

# the order of libraries matters, dependances should be accounted in left2right order
target_link_libraries(${TARGET_NAME} PRIVATE assetslib mayalib heightfieldlib corelib Qt6::Core Qt6::Gui Qt6::Widgets Qt6::Charts Qt6::PrintSupport ${GLEW_LIBRARIES} ${OPENGL_LIBRARIES})

set_target_properties(${TARGET_NAME} PROPERTIES
        #    WIN32_EXECUTABLE ON
        MACOSX_BUNDLE ON
)
