project(
  'stellgap',
  ['c', 'fortran'],
  version: '0.1',
  license: 'BSD-3',
  meson_version: '>=1.4.0',
  default_options: [
    'warning_level=2',
    'default_library=static',
    'optimization=2',
    'fortran_std=f2008',
  ],
)
fortran_compiler = meson.get_compiler('fortran')
fs = import('fs')

if fortran_compiler.get_id() == 'gcc'
  add_project_arguments('-ffree-line-length-none', language: 'fortran')
endif

# DEPENDENCIES
netcdf = dependency('netcdf', language: 'fortran')
lapack = dependency('lapack')
blas = dependency('blas')
stell = fortran_compiler.find_library(
  'stell',
  dirs: fs.expanduser('~/bin/libstell_dir/'),
)

# FORTRAN LIBRARY
lib_stgap_common_src = files(
  'src/common/constants.f90',
  'src/common/fitpack.f90',
  'src/common/fourier_lib.f90',
  'src/common/fourier_lib_convolve.f90',
  'src/common/helper.f90',
)

# lib_stgap_common = library(
#   'stgap_common',
#   lib_stgap_common_src,
#   install: false,
#   pic: true,
# )

include_dirs = include_directories(
  fs.expanduser('~/bin/libstell_dir/'),
  '/usr/include/',
)

xmetric_dependencies = [netcdf, lapack, blas, stell]
xmetric_files = files(
  'src/common/constants.f90',
  'src/metric/globals.f90',
  'src/metric/output.f90',
  'src/metric/read_input.f90',
)
xmetric = executable(
  'xmetric',
  ['src/metric/metric_element_create.f90', xmetric_files],
  include_directories: include_dirs,
  dependencies: xmetric_dependencies,
  # link_with: [lib_stgap_common],
)
