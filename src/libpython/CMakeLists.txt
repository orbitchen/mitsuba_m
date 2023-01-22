add_compile_definitions(MTS_BUILD_MODULE=MTS_MODULE_PYTHON)

add_library(mitsuba-python SHARED core.cpp render.cpp)
target_compile_options(mitsuba-python PRIVATE -fno-strict-aliasing)
set_target_properties(mitsuba-python PROPERTIES
  LIBRARY_OUTPUT_NAME mitsuba
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python${Python_VERSION}
  PREFIX "")
