# in order to propoer work, you have to set this env variable to point to your Qt install dir
# CMAKE_PREFIX_PATH
cmake_minimum_required(VERSION 3.16)
project(VectorTerrains VERSION 1.0.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_FLAGS "-mavx -mavx2 -march=native")
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets Charts PrintSupport)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOUIC ON)

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

        Source/libs/analyticfield.cpp
        Source/libs/analyticpalette.cpp
        Source/libs/analyticscalarfield.cpp
        Source/libs/array.cpp
        Source/libs/array2.cpp
        Source/libs/axis.cpp
        Source/libs/box.cpp
        Source/libs/box2.cpp
        Source/libs/camera.cpp
        Source/libs/capsule.cpp
        Source/libs/chart.cpp
        Source/libs/circle.cpp
        Source/libs/circle2.cpp
        Source/libs/color-hsl.cpp
        Source/libs/color.cpp
        Source/libs/cone.cpp
        Source/libs/convex2.cpp
        Source/libs/cpu.cpp
        Source/libs/cubic.cpp
        Source/libs/cuboctahedron.cpp
        Source/libs/cuboid.cpp
        Source/libs/curvecubic.cpp
        Source/libs/curvecubic2.cpp
        Source/libs/curvequadric.cpp
        Source/libs/curvequadric2.cpp
        Source/libs/curvequadric2set.cpp
        Source/libs/curvesetcubic.cpp
        Source/libs/curvesetcubic2.cpp
        Source/libs/curvesetquadric.cpp
        Source/libs/cylinder.cpp
        Source/libs/disc2.cpp
        Source/libs/disctiles.cpp
        Source/libs/dodecahedron.cpp
        Source/libs/draw.cpp
        Source/libs/ellipse2.cpp
        Source/libs/evector.cpp
        Source/libs/evector2.cpp
        Source/libs/evectorfloat.cpp
        Source/libs/frame.cpp
        Source/libs/frame2.cpp
        Source/libs/framescaled.cpp
        Source/libs/gpu-shader.cpp
        Source/libs/gpuanalysis.cpp
        Source/libs/gpuheightfield.cpp
        Source/libs/gputerrainsimu.cpp
        Source/libs/gputerrainsimu_test.cpp
        Source/libs/hashtable.cpp
        Source/libs/heightfield-access.cpp
        Source/libs/heightfield-carve.cpp
        Source/libs/heightfield-flow.cpp
        Source/libs/heightfield-mesh.cpp
        Source/libs/heightfield.cpp
        Source/libs/heightfieldshader.cpp
        Source/libs/hemisphere.cpp
        Source/libs/hexagon2.cpp
        Source/libs/histogram.cpp
        Source/libs/ia.cpp
        Source/libs/icosahedron.cpp
        Source/libs/icosidodecahedron.cpp
        Source/libs/integerfield2.cpp
        Source/libs/line2.cpp
        Source/libs/linear.cpp
        Source/libs/lines2.cpp
        Source/libs/mathematics.cpp
        Source/libs/matrix.cpp
        Source/libs/matrix2.cpp
        Source/libs/matrix4.cpp
        Source/libs/maya-cameraset.cpp
        Source/libs/maya-geometry.cpp
        Source/libs/maya-geometryset.cpp
        Source/libs/maya-glwidget.cpp
        Source/libs/maya-instancevbo.cpp
        Source/libs/maya-instancevboall.cpp
        Source/libs/maya-instancevboset.cpp
        Source/libs/maya-material.cpp
        Source/libs/maya-plane.cpp
        Source/libs/maya-renderer.cpp
        Source/libs/maya-shader.cpp
        Source/libs/maya-statistics.cpp
        Source/libs/maya-widget-initialize.cpp
        Source/libs/mesh-widget.cpp
        Source/libs/mesh.cpp
        Source/libs/mesh2.cpp
        Source/libs/meshcolor.cpp
        Source/libs/noise-simplex.cpp
        Source/libs/noise.cpp
        Source/libs/nonic.cpp
        Source/libs/octahedron.cpp
        Source/libs/octic.cpp
        Source/libs/octogon2.cpp
        Source/libs/palette.cpp
        Source/libs/particle2.cpp
        Source/libs/pentagon2.cpp
        Source/libs/plane.cpp
        Source/libs/pointcurve.cpp
        Source/libs/pointcurve2.cpp
        Source/libs/polygon.cpp
        Source/libs/polygon2.cpp
        Source/libs/polygons2.cpp
        Source/libs/polynomial-solver.cpp
        Source/libs/polynomial.cpp
        Source/libs/pyramid.cpp
        Source/libs/quad.cpp
        Source/libs/quad2.cpp
        Source/libs/quadric.cpp
        Source/libs/quadricsurface.cpp
        Source/libs/quartic.cpp
        Source/libs/quaternion.cpp
        Source/libs/quintic.cpp
        Source/libs/random.cpp
        Source/libs/randomfast.cpp
        Source/libs/ray.cpp
        Source/libs/rectangle.cpp
        Source/libs/sampling2.cpp
        Source/libs/scalarfield-mesh.cpp
        Source/libs/scalarfield-wyvill.cpp
        Source/libs/scalarfield.cpp
        Source/libs/scalarfield2-curvature.cpp
        Source/libs/scalarfield2-morph.cpp
        Source/libs/scalarfield2-smooth.cpp
        Source/libs/scalarfield2.cpp
        Source/libs/segment.cpp
        Source/libs/segment2.cpp
        Source/libs/septic.cpp
        Source/libs/sextic.cpp
        Source/libs/shuffle.cpp
        Source/libs/smoothellipse2.cpp
        Source/libs/sphere.cpp
        Source/libs/sphereset.cpp
        Source/libs/spheretracingwidget.cpp
        Source/libs/system.cpp
        Source/libs/terrain-raytracingwidget.cpp
        Source/libs/tetrahedra.cpp
        Source/libs/torus.cpp
        Source/libs/triangle.cpp
        Source/libs/triangle2.cpp
        Source/libs/turbulence.cpp
        Source/libs/vectorfield.cpp
        Source/libs/vectorfield2.cpp
        Source/libs/vectorset.cpp
        Source/libs/vectorset2.cpp
        Source/libs/voxel.cpp
)

SET(UI_FILES
        UI/main.ui
)

qt_wrap_ui(SOURCE_FILES UI/main.ui)
qt_wrap_cpp(SOURCE_FILES Include/qte.h)
qt_wrap_cpp(SOURCE_FILES Include/GaussianTerrainRaytracingWidget.h Include/libs/maya.h Include/libs/realtime.h)

qt_add_executable(${PROJECT_NAME}
        ${SOURCE_FILES}
        ${UI_FILES}
)

target_compile_definitions(${PROJECT_NAME} PUBLIC SOLUTION_DIR="${CMAKE_CURRENT_SOURCE_DIR}/")
target_include_directories(${PROJECT_NAME} PRIVATE Include/)
target_include_directories(${PROJECT_NAME} PRIVATE ${GLEW_INCLUDE_DIRS})

target_link_directories(${PROJECT_NAME} BEFORE PUBLIC $ENV{GLEW_INCLUDE_DIRS})

target_link_libraries(${PROJECT_NAME} PRIVATE Qt6::Core Qt6::Gui Qt6::Widgets Qt6::Charts Qt6::PrintSupport ${GLEW_LIBRARIES} ${OPENGL_LIBRARIES})

set_target_properties(${PROJECT_NAME} PROPERTIES
        #    WIN32_EXECUTABLE ON
        MACOSX_BUNDLE ON
)
