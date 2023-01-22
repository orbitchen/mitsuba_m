cmake_minimum_required(VERSION 3.16)

project(Mitsuba
  VERSION 0.6.0
  DESCRIPTION
    "Personal Compilable Version of the Mitsuba Renderer"
  LANGUAGES
    CXX C)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build Mitsuba in ${CMAKE_BUILD_TYPE} mode")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS "-msse2 -fopenmp -fPIC")

include_directories(include)
link_libraries(stdc++ dl m pthread gomp)
add_compile_definitions(GLEW_MX MTS_DEBUG SINGLE_PRECISION SPECTRUM_SAMPLES=3 MTS_SSE MTS_HAS_COHERENT_RT)

# Check necessary dependencies
include(FindPkgConfig)

find_package(XercesC REQUIRED)
link_libraries(${XercesC_LIBRARIES})

find_package(Boost REQUIRED COMPONENTS system filesystem thread python)
link_libraries(Boost::system Boost::filesystem Boost::thread Boost::python)

find_package(Eigen3 REQUIRED CONFIG)
include_directories(${EIGEN3_INCLUDE_DIRS})

find_package(OpenMP REQUIRED)

set(OpenGL_GL_PREFERENCE LEGACY)
find_package(OpenGL REQUIRED)
link_libraries(OpenGL::GL OpenGL::GLU)

pkg_check_modules(GLEWMX REQUIRED glewmx)
link_libraries(${GLEWMX_LIBRARIES})

pkg_check_modules(Xxf86vm REQUIRED xxf86vm)
link_libraries(${Xxf86vm_LIBRARIES})

find_package(X11 REQUIRED)
link_libraries(${X11_LIBRARIES})

# Check optional Dependencies
find_package(PNG)
if(NOT PNG_FOUND)
  message(WARNING "PNG library not found -- PNG I/O not available")
else()
  add_compile_definitions(MTS_HAS_LIBPNG=1)
  link_libraries(${PNG_LIBRARIES})
endif()

find_package(JPEG)
if(NOT JPEG_FOUND)
  message(WARNING "JPEG library not found -- JPEG I/O not available")
else()
  add_compile_definitions(MTS_HAS_LIBJPEG=1)
  link_libraries(${JPEG_LIBRARIES})
endif()

pkg_check_modules(OpenEXR OpenEXR)
if(NOT OpenEXR_FOUND)
  message(WARNING "OpenEXR library not found -- OpenEXR I/O not available")
else()
  add_compile_definitions(MTS_HAS_OPENEXR=1)
  include_directories(${OpenEXR_INCLUDE_DIRS})
  link_libraries(${OpenEXR_LIBRARIES})
endif()

find_package(COLLADA_DOM COMPONENTS 1.4.1 CONFIG)
if(NOT COLLADA_DOM_FOUND)
  message(WARNING "Collada library not found -- Collada importer is disabled")
else()
  add_compile_definitions(MTS_HAS_COLLADA=1)
  include_directories(${COLLADA_DOM_INCLUDE_DIRS})
  link_libraries(${COLLADA_DOM_LIBRARIES})
endif()

pkg_check_modules(FFTW3 fftw3)
if(NOT FFTW3_FOUND)
  message(WARNING "FFTW3 library not found -- Fast image convolution not available")
else()
  add_compile_definitions(MTS_HAS_FFTW=1)
  link_libraries(fftw3_threads fftw3)
endif()

find_package(Qt5 COMPONENTS Gui Widgets Core OpenGL Xml XmlPatterns Network)
if(NOT Qt5_FOUND)
  message(WARNING "Qt5 library not found -- Mitsuba GUI not available")
else()
  link_libraries(Qt5::Gui Qt5::Widgets Qt5::Core Qt5::OpenGL Qt5::Xml Qt5::XmlPatterns Qt5::Network)
endif()

find_package(Python COMPONENTS Interpreter Development)
if(NOT Python_FOUND)
  message(WARNING "Python interpreter or library not found -- Python API not available")
else()
  include_directories(${Python_INCLUDE_DIRS})
  link_libraries(${Python_LIBRARIES})
endif()

# Build the source code
add_subdirectory(src)

# Copy the data files needed at runtime
file(COPY data/schema DESTINATION ${CMAKE_BINARY_DIR}/data REGEX ".*\\.xsd$|.*\\.xsl$")
file(COPY data/microfacet DESTINATION ${CMAKE_BINARY_DIR}/data REGEX ".*\\.dat$")
file(COPY data/ior DESTINATION ${CMAKE_BINARY_DIR}/data REGEX ".*\\.spd$")
