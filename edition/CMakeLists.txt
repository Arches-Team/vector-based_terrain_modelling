# Tested on Windows - in order to proper work, you have to set this env variable to point to your Qt install dir
# CMAKE_PREFIX_PATH=Qt/versionx.x/compiler/lib/cmake
cmake_minimum_required(VERSION 3.16)
project(VectorTerrains VERSION 1.0.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_FLAGS "-mavx -mavx2 -march=native")
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets Charts PrintSupport)

# OpenGL & GLEW library
find_package(OpenGL)
find_package(GLEW)
if (NOT GLEW_FOUND)
    message("Adding GLEW manually for Windows...")
    set(PLATFORM x64)
    set(GLEW_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include/libs/GL/include)
    set(GLEW_LIBRARIES ${CMAKE_CURRENT_SOURCE_DIR}/lib/$(PLATFORM)/glew32.lib)
    set(GLEW_FOUND true)
endif ()

set(SOURCE_FILES
        src/GaussianTerrainRaytracingWidget.cpp
        src/graph.cpp
        src/Kernels.cpp
        src/main.cpp
        src/MayaSimpleRendererColors.cpp
        src/qte-main-window.cpp
        src/utils.cpp

        src/Tools/ToolBrush.cpp
        src/Tools/ToolEdit.cpp
        src/Tools/ToolGraph.cpp

        src/Kernels/DetailsKernel.cpp
        src/Kernels/GaussianKernel.cpp
        src/Kernels/Kernel.cpp

        src/libs/analyticpalette.cpp
        src/libs/analyticscalarfield.cpp
        src/libs/array.cpp
        src/libs/array2.cpp
        src/libs/axis.cpp
        src/libs/box.cpp
        src/libs/box2.cpp
        src/libs/camera.cpp
        src/libs/capsule.cpp
        src/libs/chart.cpp
        src/libs/circle.cpp
        src/libs/circle2.cpp
        src/libs/color-hsl.cpp
        src/libs/color.cpp
        src/libs/cone.cpp
        src/libs/convex2.cpp
        src/libs/cpu.cpp
        src/libs/cubic.cpp
        src/libs/cuboctahedron.cpp
        src/libs/cuboid.cpp
        src/libs/curvecubic.cpp
        src/libs/curvecubic2.cpp
        src/libs/curvequadric.cpp
        src/libs/curvequadric2.cpp
        src/libs/curvequadric2set.cpp
        src/libs/curvesetcubic.cpp
        src/libs/curvesetcubic2.cpp
        src/libs/curvesetquadric.cpp
        src/libs/cylinder.cpp
        src/libs/disc2.cpp
        src/libs/disctiles.cpp
        src/libs/dodecahedron.cpp
        src/libs/draw.cpp
        src/libs/ellipse2.cpp
        src/libs/evector.cpp
        src/libs/evector2.cpp
        src/libs/evectorfloat.cpp
        src/libs/frame.cpp
        src/libs/frame2.cpp
        src/libs/framescaled.cpp
        src/libs/gpu-shader.cpp
        src/libs/gpuanalysis.cpp
        src/libs/gpuheightfield.cpp
        src/libs/gputerrainsimu.cpp
        src/libs/gputerrainsimu_test.cpp
        src/libs/hashtable.cpp
        src/libs/heightfield-access.cpp
        src/libs/heightfield-carve.cpp
        src/libs/heightfield-flow.cpp
        src/libs/heightfield-mesh.cpp
        src/libs/heightfield.cpp
        src/libs/heightfieldshader.cpp
        src/libs/hemisphere.cpp
        src/libs/hexagon2.cpp
        src/libs/histogram.cpp
        src/libs/ia.cpp
        src/libs/icosahedron.cpp
        src/libs/icosidodecahedron.cpp
        src/libs/integerfield2.cpp
        src/libs/line2.cpp
        src/libs/linear.cpp
        src/libs/lines2.cpp
        src/libs/mathematics.cpp
        src/libs/matrix.cpp
        src/libs/matrix2.cpp
        src/libs/matrix4.cpp
        src/libs/maya-cameraset.cpp
        src/libs/maya-geometry.cpp
        src/libs/maya-geometryset.cpp
        src/libs/maya-glwidget.cpp
        src/libs/maya-instancevbo.cpp
        src/libs/maya-instancevboall.cpp
        src/libs/maya-instancevboset.cpp
        src/libs/maya-material.cpp
        src/libs/maya-plane.cpp
        src/libs/maya-renderer.cpp
        src/libs/maya-shader.cpp
        src/libs/maya-statistics.cpp
        src/libs/maya-widget-initialize.cpp
        src/libs/mesh-widget.cpp
        src/libs/mesh.cpp
        src/libs/mesh2.cpp
        src/libs/meshcolor.cpp
        src/libs/noise-simplex.cpp
        src/libs/noise.cpp
        src/libs/nonic.cpp
        src/libs/octahedron.cpp
        src/libs/octic.cpp
        src/libs/octogon2.cpp
        src/libs/palette.cpp
        src/libs/particle2.cpp
        src/libs/pentagon2.cpp
        src/libs/plane.cpp
        src/libs/pointcurve.cpp
        src/libs/pointcurve2.cpp
        src/libs/polygon.cpp
        src/libs/polygon2.cpp
        src/libs/polygons2.cpp
        src/libs/polynomial-solver.cpp
        src/libs/polynomial.cpp
        src/libs/pyramid.cpp
        src/libs/quad.cpp
        src/libs/quad2.cpp
        src/libs/quadric.cpp
        src/libs/quadricsurface.cpp
        src/libs/quartic.cpp
        src/libs/quaternion.cpp
        src/libs/quintic.cpp
        src/libs/random.cpp
        src/libs/randomfast.cpp
        src/libs/ray.cpp
        src/libs/rectangle.cpp
        src/libs/sampling2.cpp
        src/libs/scalarfield-mesh.cpp
        src/libs/scalarfield-wyvill.cpp
        src/libs/scalarfield.cpp
        src/libs/scalarfield2-curvature.cpp
        src/libs/scalarfield2-morph.cpp
        src/libs/scalarfield2-smooth.cpp
        src/libs/scalarfield2.cpp
        src/libs/segment.cpp
        src/libs/segment2.cpp
        src/libs/septic.cpp
        src/libs/sextic.cpp
        src/libs/shuffle.cpp
        src/libs/smoothellipse2.cpp
        src/libs/sphere.cpp
        src/libs/sphereset.cpp
        src/libs/spheretracingwidget.cpp
        src/libs/system.cpp
        src/libs/terrain-raytracingwidget.cpp
        src/libs/tetrahedra.cpp
        src/libs/torus.cpp
        src/libs/triangle.cpp
        src/libs/triangle2.cpp
        src/libs/turbulence.cpp
        src/libs/vectorfield.cpp
        src/libs/vectorfield2.cpp
        src/libs/vectorset.cpp
        src/libs/vectorset2.cpp
        src/libs/voxel.cpp
)

SET(UI_FILES
        UI/main.ui
)

qt_wrap_ui(SOURCE_FILES UI/main.ui)
qt_wrap_cpp(SOURCE_FILES include/qte.h)
qt_wrap_cpp(SOURCE_FILES include/GaussianTerrainRaytracingWidget.h include/libs/maya.h include/libs/realtime.h)

qt_add_executable(${PROJECT_NAME}
        ${SOURCE_FILES}
        ${UI_FILES}
)

target_compile_definitions(${PROJECT_NAME} PUBLIC SOLUTION_DIR="${CMAKE_CURRENT_SOURCE_DIR}/")
target_include_directories(${PROJECT_NAME} PRIVATE include/)
target_include_directories(${PROJECT_NAME} PRIVATE ${GLEW_INCLUDE_DIRS})

target_link_directories(${PROJECT_NAME} BEFORE PUBLIC $ENV{GLEW_INCLUDE_DIRS})

target_link_libraries(${PROJECT_NAME} PRIVATE Qt6::Core Qt6::Gui Qt6::Widgets Qt6::Charts Qt6::PrintSupport ${GLEW_LIBRARIES} ${OPENGL_LIBRARIES})

set_target_properties(${PROJECT_NAME} PROPERTIES
        #    WIN32_EXECUTABLE ON
        MACOSX_BUNDLE ON
)
