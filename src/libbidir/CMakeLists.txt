add_compile_definitions(MTS_BUILD_MODULE=MTS_MODULE_BIDIR)

add_mts_library(mitsuba-bidir
  common.cpp rsampler.cpp vertex.cpp edge.cpp 
  path.cpp verification.cpp util.cpp pathsampler.cpp 
  mut_bidir.cpp mut_lens.cpp mut_caustic.cpp 
  mut_mchain.cpp manifold.cpp mut_manifold.cpp)
