add_compile_definitions(COLLADA_DOM_2_4)

add_mts_library(mitsuba-converter SHARED collada.cpp obj.cpp converter.cpp)
