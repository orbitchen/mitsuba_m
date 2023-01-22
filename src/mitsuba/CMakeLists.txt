add_mts_executable(mtsutil mtsutil.cpp)
add_mts_executable(mtssrv mtssrv.cpp)
add_mts_executable(mitsuba mitsuba.cpp)

if(COLLADA_DOM_FOUND)
  add_mts_executable(mtsimport mtsimport.cpp)
endif()
