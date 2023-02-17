add_library(oscalib plinklite.c bodfile.c)
target_include_directories(oscalib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
