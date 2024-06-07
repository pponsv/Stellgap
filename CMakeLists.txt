cmake_minimum_required(VERSION 3.15)

enable_language(Fortran)

project(Stellgap)

get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

# Debug:

# set(CMAKE_VERBOSE_MAKEFILE ON)

# Options
set(SERIAL "SERIAL")

# Paths
set(LIBSTELL_PATH "$ENV{HOME}/bin/libstell_dir/")
set(NETCDF_PATH "/usr/include")

# Libraries
add_library(libstell SHARED IMPORTED)
set_target_properties(libstell PROPERTIES IMPORTED_LOCATION ${LIBSTELL_PATH}/libstell.so)

# XMETRIC
set(XMETRIC_SRC "src_old/metric_element_create_ver8.46.f90" "src_old/constants.f90")
add_executable(xmetric ${XMETRIC_SRC})
target_include_directories(xmetric PRIVATE ${LIBSTELL_PATH})
target_link_directories(xmetric PRIVATE ${LIBSTELL_PATH})
target_link_libraries(xmetric libstell netcdf lapack blas)

# XSTGAP
set(XSTGAP_SRC "src_old/stellgap_ver5.f90")
set_source_files_properties(
  ${XSTGAP_SRC}
  PROPERTIES Fortran_PREPROCESS ON
)
set(MODULES_SRC "src_old/constants.f90" "src_old/fourier_lib.f90" "src_old/Fourier_lib_convolve.f90" "src_old/post_process.f90" "src_old/fitpack.f90" "src_old/read_input.f90" "src_old/output.f90")
add_executable(xstgap ${XSTGAP_SRC} ${MODULES_SRC})
target_link_libraries(xstgap lapack blas m pthread)
target_compile_definitions(xstgap PRIVATE -D${SERIAL})
target_compile_options(xstgap PRIVATE -cpp -P -C -O2 -ffree-form)

# XSTGAP_SND
set(XSTGAP_SND_SRC "src_old/stellgap_soundwave_lagrng_ver7.f90")
set_source_files_properties(
  ${XSTGAP_SND_SRC}
  PROPERTIES Fortran_PREPROCESS ON
)
set(MODULES_SRC "src_old/constants.f90" "src_old/fourier_lib.f90" "src_old/Fourier_lib_convolve.f90" "src_old/post_process.f90" "src_old/fitpack.f90" "src_old/read_input.f90" "src_old/output.f90")
add_executable(xstgap_snd ${XSTGAP_SND_SRC} ${MODULES_SRC})
target_link_libraries(xstgap_snd lapack blas m pthread)
target_compile_definitions(xstgap_snd PRIVATE -D${SERIAL})
target_compile_options(xstgap_snd PRIVATE -cpp -P -C -O2 -ffree-form)

# XSTGAP_NEW
set(XSTGAP_NEW_SRC "src/stellgap/stellgap_ver5_new.f90")
set_source_files_properties(
  ${XSTGAP_NEW_SRC}
  PROPERTIES Fortran_PREPROCESS ON
)
set(MODULES_COMMON "src/common/constants.f90" "src/common/fourier_lib.f90" "src/common/fourier_lib_convolve.f90" "src/common/fitpack.f90" "src/common/helper.f90")
set(MODULES_STELLGAP "src/stellgap/post_process.f90" "src/stellgap/read_input.f90" "src/stellgap/output.f90" "src/stellgap/globals.f90" "src/stellgap/eig_solver.f90")
add_executable(xstgap_new ${XSTGAP_NEW_SRC} ${MODULES_COMMON} ${MODULES_STELLGAP})
target_include_directories(xstgap_new PRIVATE ${NETCDF_PATH})
target_link_libraries(xstgap_new lapack blas m pthread netcdff)
target_compile_definitions(xstgap_new PRIVATE -D${SERIAL})
target_compile_options(xstgap_new PRIVATE -cpp -P -C -O3 -ffree-form)

# XMETRIC_NEW
set(MODULES_XMETRIC_COMMON "src/common/constants.f90")
set(MODULES_XMETRIC "src/metric/read_input.f90" "src/metric/output.f90" "src/metric/metric_globals.f90")
set(XMETRIC_SRC "src/metric/metric_element_create.f90" ${MODULES_XMETRIC_COMMON} ${MODULES_XMETRIC})
add_executable(xmetric_new ${XMETRIC_SRC})
target_include_directories(xmetric_new PRIVATE ${LIBSTELL_PATH})
target_link_directories(xmetric_new PRIVATE ${LIBSTELL_PATH})
target_link_libraries(xmetric_new libstell netcdf lapack blas)
