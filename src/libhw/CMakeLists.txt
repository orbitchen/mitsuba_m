add_compile_definitions(MTS_BUILD_MODULE=MTS_MODULE_HW)

add_mts_library(mitsuba-hw
  session.cpp device.cpp gputexture.cpp gpugeometry.cpp
  gpuprogram.cpp renderer.cpp glrenderer.cpp glprogram.cpp
  glgeometry.cpp gltexture.cpp gpusync.cpp glsync.cpp
  vpl.cpp font.cpp viewer.cpp basicshader.cpp shadow.cpp)

if(${CMAKE_SYSTEM_NAME} MATCHES Linux)
  target_sources(mitsuba-hw PRIVATE
    x11session.cpp x11device.cpp glxdevice.cpp glxrenderer.cpp)
endif()
