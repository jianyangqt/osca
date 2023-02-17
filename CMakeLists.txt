#osca CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(osca C CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)


#set project source directory
set(PRJ_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
set(CUSTORM_INCLUDE_SERACH_PATH "" CACHE PATH "")
set(CUSTORM_LIB_SEARCH_PATH "" CACHE PATH "")
set(BUILD_STATIC FALSE CACHE BOOL "")


#information of CPU architecture
message(STATUS "CPU_ARCH: " ${CMAKE_HOST_SYSTEM_PROCESSOR})
message("")

#information of operation system
message(STATUS ${CMAKE_HOST_SYSTEM})
if (UNIX)
    set(system "UNIX-like")
elseif(APPLE)
    set(system "macOS")
elseif(WIN32)
    set(system "windows")
else()
    message(FATAL_ERROR "Operation system not recognized, only support unix apple and windows.")
endif()
message(STATUS "OS: " ${system})
message("")

#information of compilation toolchain
message(STATUS "CMAKE_C_COMPILER:  " ${CMAKE_C_COMPILER} " " ${CMAKE_C_COMPILER_ID} ":" ${CMAKE_C_COMPILER_VERSION} " " ${CMAKE_C_STANDARD})
message(STATUS "CMAKE_CXX_COMPILER: " ${CMAKE_CXX_COMPILER} " " ${CMAKE_CXX_COMPILER_ID} ":" ${CMAKE_CXX_COMPILER_VERSION} " " ${CMAKE_CXX_STANDARD})
message(STATUS "CMAKE_AR: " ${CMAKE_AR})
message(STATUS "CMAKE_RANLIB: " ${CMAKE_RANLIB})
message(STATUS "CMAKE_COMMAND: " ${CMAKE_COMMAND})
message(STATUS "CMAKE_VERSION: " ${CMAKE_VERSION})
message(STATUS "CMAKE_GENERATOR: " ${CMAKE_GENERATOR})
message(STATUS "CMAKE_BUILD_TOOL: " ${CMAKE_BUILD_TOOL})

message(STATUS "CMAKE_SHARED_LIBRARY_PREFIX: " ${CMAKE_SHARED_LIBRARY_PREFIX})
message(STATUS "CMAKE_SHARED_LIBRARY_SUFFIX: " ${CMAKE_SHARED_LIBRARY_SUFFIX})
message(STATUS "CMAKE_STATIC_LIBRARY_PREFIX: " ${CMAKE_STATIC_LIBRARY_PREFIX})
message(STATUS "CMAKE_STATIC_LIBRARY_SUFFIX: " ${CMAKE_STATIC_LIBRARY_SUFFIX})

message(STATUS "CMAKE_FIND_LIBRARY_PREFIXES: " ${CMAKE_FIND_LIBRARY_PREFIXES})
#set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
message(STATUS "CMAKE_FIND_LIBRARY_SUFFIXES: " ${CMAKE_FIND_LIBRARY_SUFFIXES})

message(STATUS "CMAKE_PROJECT_NAME: " ${CMAKE_PROJECT_NAME})
message(STATUS "CMAKE_PROJECT_VERSION: " ${CMAKE_PROJECT_VERSION})
message(STATUS "CMAKE_BINARY_DIR: " ${CMAKE_BINARY_DIR})
message(STATUS "CMAKE_SOURCE_DIR: " ${CMAKE_SOURCE_DIR})
message("")


#detect zlib
set(zlib_path FALSE CACHE PATH "zlib package location")
set(zlib_include_path FALSE CACHE PATH "zlib head file location")
set(zlib_lib_path FALSE CACHE PATH "zlib library location")
if (zlib_path)
    set(zlib_include_dir ${zlib_path}/include)
    set(zlib_lib_dir ${zlib_path}/lib)
elseif (zlib_include_path AND zlib_lib_path)
    set(zlib_include_dir ${zlib_include_path})
    set(zlib_lib_dir ${zlib_lib_path})
else()
    find_path(zlib_include_dir zlib.h ${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES} ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
    find_library(zlib_lib_dir z ${CMAKE_C_IMPLICIT_LINK_DIRECTORIES} ${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES})
    if (zlib_include_dir AND zlib_lib_dir)
        get_filename_component(zlib_lib_dir ${zlib_lib_dir} DIRECTORY)
        message(STATUS "zlib founded.")
    else()
        message(FATAL_ERROR "Can not find zlib lib, please specific them by -Dzlib_path or -Dzlib_include_path and -Dzlib_lib_path")
    endif()
endif()


#detect gsl
set(gsl_path FALSE CACHE PATH "gsl package directory")
set(gsl_include_path FALSE CACHE PATH "gsl include direcoty")
set(gsl_lib_path FALSE CACHE PATH "gsl lib directory")
if (gsl_path)
    set(gsl_include_dir ${gsl_path}/include)
    set(gsl_lib_dir ${gsl_path}/lib)
    message(STATUS "Set gsl head files: " ${gsl_include_dir})
    message(STATUS "Set gsl lib files: " ${gsl_lib_dir})
elseif (gsl_include_path AND gsl_lib_path)
    set(gsl_include_dir ${gsl_include_path})
    set(gsl_lib_dir ${gsl_lib_path})
    message(STATUS "Set gsl head files: " ${gsl_include_dir})
    message(STATUS "Set gsl lib files: " ${gsl_lib_dir})
else()
    find_path(gsl_include_dir gsl/gsl_sf_bessel.h  ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES} ${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES})
    find_library(gsl_lib_dir gsl ${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES} ${CMAKE_C_IMPLICIT_LINK_DIRECTORIES})
    if (gsl_include_dir AND gsl_lib_dir)
        get_filename_component(gsl_lib_dir ${gsl_lib_dir} DIRECTORY)
        message(STATUS "Find gsl head files: " ${gsl_include_dir})
        message(STATUS "Find gls lib files: " ${gsl_lib_dir})
    else ()
        message(SEND_ERROR "Not find gsl head files and libs.")
        message(SEND_ERROR "please specific them by gsl_dir or gsl_include_dir and gls_lib_dir")
        return()
    endif()
endif()


#detect Eigen
set(eigen_path FALSE CACHE PATH "path where Eigen located.")
#set(eigen_path FALSE)
if (eigen_path)
    message(STATUS "You set eigen path to: " ${eigen_path})
    set(eigen_dir ${eigen_path}/eigen3)
else()
    find_path(eigen_dir "eigen3/Eigen" ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES} ${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES})
    if (eigen_dir)
        set(eigen_dir ${eigen_dir}/eigen3)
        message(STATUS "Find Eigen: " ${eigen_dir})
    else()
        message(FATAL_ERROR "can not find eigen, please using -Deigen_path=/where/your/eigen/located to set eigen path")
    endif()
endif()


#detect intel MKL
set(mkl_path FALSE CACHE PATH "intel MKL package location")
set(mkl_include_path CACHE PATH "")
set(mkl_lib_path CACHE PATH "")
if (mkl_path)
    message(STATUS "You set mkl path to: " ${mkl_path})
    set(mkl_include_dir ${mkl_path}/include)
    set(mkl_lib_dir ${mkl_path}/lib/intel64)
elseif (mkl_include_path AND mkl_lib_path)
    message(STATUS "You set mkl_includt_path to: " ${mkl_include_path})
    message(STATUS "You set mkl_lib_path to: " ${mkl_lib_path})
    set(mkl_include_dir ${mkl_include_path})
    set(mkl_lib_dir ${mkl_lib_path}/intel64)
else ()
    find_path(mkl_include_dir mkl.h ${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES} ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
    find_library(mkl_lib_dir mkl_core ${CMAKE_C_IMPLICIT_LINK_DIRECTORIES} ${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES})
    if (mkl_include_dir AND mkl_lib_dir)
        get_filename_component(mkl_lib_dir ${mkl_lib_dir} DIRECTORY)
        message("Find mkl head file and ")
    else()
        message(FATAL_ERROR "Can not find mkl head files and libs, please specific them by -Dmkl_path or -Dmkl_include_path and -Dmkl_lib_path to specific them.")
    endif()
endif()


#detech Rmath
set(rmath_path FALSE CACHE PATH "")
set(rmath_include_path CACHE PATH "")
set(rmath_lib_path CACHE PATH "")
if (rmath_path)
    message(STATUS "You set rmath_path to: " ${rmath_path})
    set(rmath_include_dir ${rmath_path}/include)
    set(rmath_lib_dir ${rmath_path}/lib)
elseif (rmath_include_path AND rmath_lib_path)
    message(STATUS "You set rmath_include_path to: " ${rmath_include_path})
    message(STATUS "You set rmath_lib_path to: " ${rmath_lib_path})
    set(rmath_include_dir ${rmath_include_path})
    set(rmath_lib_dir ${rmath_lib_path})
else()
    find_path(rmath_include_dir Rmath.h ${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES} ${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES})
    find_library(rmath_lib_dir Rmath ${CMAKE_C_IMPLICIT_LINK_DIRECTORIES} ${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES})
    if (rmath_include_dir AND rmath_lib_dir)
        get_filename_component(rmath_lib_dir ${rmath_lib_dir} DIRECTORY)
        message(STATUS "Find Rmath head files at: " ${rmath_include_dir})
        message(STATUS "Find Rmath lib files at: " ${rmath_lib_dir})
    else()
        message(FATAL_ERROR "can not find Rmath head files and libs, pliease specific them by -Drmath_path or -Drmath_include_path and -Drmath_lib_path")
    endif()
endif()


#build
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
#set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
add_executable(osca ${PRJ_SOURCE_DIR}/l0_com.cpp ${PRJ_SOURCE_DIR}/l0_io.cpp 
    ${PRJ_SOURCE_DIR}/l0_mem.cpp ${PRJ_SOURCE_DIR}/l0_stat.cpp 
    ${PRJ_SOURCE_DIR}/l1_op_geno.cpp ${PRJ_SOURCE_DIR}/l1_stat.cpp
    ${PRJ_SOURCE_DIR}/l2_besd.cpp  ${PRJ_SOURCE_DIR}/l2_bfile.cpp 
    ${PRJ_SOURCE_DIR}/l2_efile.cpp  ${PRJ_SOURCE_DIR}/l2_enet.cpp
    ${PRJ_SOURCE_DIR}/l2_reml.cpp  ${PRJ_SOURCE_DIR}/l3_efile.cpp  
    ${PRJ_SOURCE_DIR}/l3_ewas.cpp  ${PRJ_SOURCE_DIR}/l3_glmnet.cpp  
    ${PRJ_SOURCE_DIR}/l3_gwas.cpp  ${PRJ_SOURCE_DIR}/l3_permutation.cpp  
    ${PRJ_SOURCE_DIR}/l3_smr.cpp  ${PRJ_SOURCE_DIR}/l3_vqtl.cpp  
    ${PRJ_SOURCE_DIR}/l4_osc.cpp  ${PRJ_SOURCE_DIR}/cis_learn_beta.cpp  
    ${PRJ_SOURCE_DIR}/dcdflib.cpp ${PRJ_SOURCE_DIR}/Module_vqtl_drm_svlm.c
    )
add_subdirectory(lib)

if (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "x86_64")
    if (system STREQUAL "UNIX-like")
        if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
            target_compile_options(osca PRIVATE -fopenmp)
            target_include_directories(osca PRIVATE ${zlib_include_dir} PRIVATE ${gsl_include_dir} PRIVATE ${eigen_dir} PRIVATE  ${mkl_include_dir} PRIVATE ${rmath_lib_dir})
            target_link_directories(osca PRIVATE ${zlib_lib_dir} PRIVATE ${gsl_lib_dir} PRIVATE ${mkl_lib_dir} PRIVATE  ${rmath_lib_dir})
            if (NOT BUILD_STATIC)
                target_link_libraries(osca z -Wl,--start-group mkl_core mkl_intel_lp64 mkl_gnu_thread -Wl,--end-group gomp pthread Rmath gsl gslcblas m oscalib)
            else()
                message(WARNING "You now try to link libs staticlly, You need make sure all static lib was installed.")
                target_link_libraries(osca -static z -Wl,--start-group mkl_core mkl_intel_lp64 mkl_gnu_thread -Wl,--end-group gomp pthread Rmath gsl gslcblas m dl)
            endif()
        elseif (CMAKE_C_COMPILER_ID STREQUAL "Clang")
            target_include_directories(osca PRIVATE ${zlib_include_dir} PRIVATE ${gsl_include_dir} PRIVATE ${eigen_dir} PRIVATE  ${mkl_include_dir} PRIVATE ${rmath_lib_dir})
            target_link_directories(osca PRIVATE ${zlib_lib_dir} PRIVATE ${gsl_lib_dir} PRIVATE ${mkl_lib_dir} PRIVATE  ${rmath_lib_dir})
            if (NOT BUILD_STATIC)
                target_link_libraries(osca z -Wl,--start-group mkl_core mkl_intel_lp64 mkl_gnu_thread -Wl,--end-group gomp pthread Rmath gsl gslcblas m)
            else()
                message(WARNING "You now try to link libs staticlly, You need make sure all static lib was installed.")
                target_link_libraries(osca -static z -Wl,--start-group mkl_core mkl_intel_lp64 mkl_gnu_thread -Wl,--end-group gomp pthread Rmath gsl gslcblas m dl)
            endif()

        endif()
    elseif (system STREQUAL "macOS")
    elseif (system STREQUAL "windows")
    else()
        message(FATAL_ERROR "OS not recgnized.")
    endif()
else()
    message(STATUS "Until now, only x86_64 architecture was tested.")
endif()


#istall
#install(TARGETS osca DESTINATION ${CMAKE_BINARY_DIR}/bin RUNTIME_DEPENDENCY_SET)
