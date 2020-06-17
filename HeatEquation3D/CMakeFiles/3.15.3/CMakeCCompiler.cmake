set(CMAKE_C_COMPILER "/apps/software/GCCcore/8.3.0/bin/cc")
set(CMAKE_C_COMPILER_ARG1 "")
set(CMAKE_C_COMPILER_ID "GNU")
set(CMAKE_C_COMPILER_VERSION "8.3.0")
set(CMAKE_C_COMPILER_VERSION_INTERNAL "")
set(CMAKE_C_COMPILER_WRAPPER "")
set(CMAKE_C_STANDARD_COMPUTED_DEFAULT "11")
set(CMAKE_C_COMPILE_FEATURES "c_std_90;c_function_prototypes;c_std_99;c_restrict;c_variadic_macros;c_std_11;c_static_assert")
set(CMAKE_C90_COMPILE_FEATURES "c_std_90;c_function_prototypes")
set(CMAKE_C99_COMPILE_FEATURES "c_std_99;c_restrict;c_variadic_macros")
set(CMAKE_C11_COMPILE_FEATURES "c_std_11;c_static_assert")

set(CMAKE_C_PLATFORM_ID "Linux")
set(CMAKE_C_SIMULATE_ID "")
set(CMAKE_C_COMPILER_FRONTEND_VARIANT "")
set(CMAKE_C_SIMULATE_VERSION "")



set(CMAKE_AR "/apps/software/binutils/2.32-GCCcore-8.3.0/bin/ar")
set(CMAKE_C_COMPILER_AR "/apps/software/GCCcore/8.3.0/bin/gcc-ar")
set(CMAKE_RANLIB "/apps/software/binutils/2.32-GCCcore-8.3.0/bin/ranlib")
set(CMAKE_C_COMPILER_RANLIB "/apps/software/GCCcore/8.3.0/bin/gcc-ranlib")
set(CMAKE_LINKER "/apps/software/binutils/2.32-GCCcore-8.3.0/bin/ld")
set(CMAKE_MT "")
set(CMAKE_COMPILER_IS_GNUCC 1)
set(CMAKE_C_COMPILER_LOADED 1)
set(CMAKE_C_COMPILER_WORKS TRUE)
set(CMAKE_C_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_C_COMPILER_ENV_VAR "CC")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_C_COMPILER_ID_RUN 1)
set(CMAKE_C_SOURCE_FILE_EXTENSIONS c;m)
set(CMAKE_C_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_C_LINKER_PREFERENCE 10)

# Save compiler ABI information.
set(CMAKE_C_SIZEOF_DATA_PTR "8")
set(CMAKE_C_COMPILER_ABI "ELF")
set(CMAKE_C_LIBRARY_ARCHITECTURE "")

if(CMAKE_C_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_C_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_C_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_C_COMPILER_ABI}")
endif()

if(CMAKE_C_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()

set(CMAKE_C_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_C_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_C_CL_SHOWINCLUDES_PREFIX}")
endif()





set(CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES "/apps/software/cURL/7.66.0-GCCcore-8.3.0/include;/apps/software/bzip2/1.0.8-GCCcore-8.3.0/include;/apps/software/ncurses/6.1-GCCcore-8.3.0/include;/apps/software/OpenMPI/3.1.4-gcccuda-2019b/include;/apps/software/hwloc/1.11.12-GCCcore-8.3.0/include;/apps/software/libpciaccess/0.14-GCCcore-8.3.0/include;/apps/software/libxml2/2.9.9-GCCcore-8.3.0/include/libxml2;/apps/software/libxml2/2.9.9-GCCcore-8.3.0/include;/apps/software/XZ/5.2.4-GCCcore-8.3.0/include;/apps/software/numactl/2.0.12-GCCcore-8.3.0/include;/apps/software/zlib/1.2.11-GCCcore-8.3.0/include;/apps/software/CUDA/10.1.243-GCC-8.3.0/nvvm/include;/apps/software/CUDA/10.1.243-GCC-8.3.0/extras/CUPTI/include;/apps/software/CUDA/10.1.243-GCC-8.3.0/include;/apps/software/imkl/2019.5.281-iimpi-2019b/mkl/include/fftw;/apps/software/imkl/2019.5.281-iimpi-2019b/mkl/include;/apps/software/impi/2018.5.288-iccifort-2019.5.281/include64;/apps/software/iccifort/2019.5.281/compilers_and_libraries_2019.5.281/linux/tbb/include;/apps/software/iccifort/2019.5.281/compilers_and_libraries_2019.5.281/linux/mkl/include;/apps/software/iccifort/2019.5.281/compilers_and_libraries_2019.5.281/linux/ipp/include;/apps/software/iccifort/2019.5.281/compilers_and_libraries_2019.5.281/linux/daal/include;/apps/software/binutils/2.32-GCCcore-8.3.0/include;/apps/software/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include;/apps/software/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include-fixed;/usr/include")
set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "gcc;gcc_s;c;gcc;gcc_s")
set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "/apps/software/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0;/apps/software/GCCcore/8.3.0/lib/gcc;/apps/software/GCCcore/8.3.0/lib64;/lib64;/usr/lib64;/apps/software/cURL/7.66.0-GCCcore-8.3.0/lib;/apps/software/bzip2/1.0.8-GCCcore-8.3.0/lib;/apps/software/ncurses/6.1-GCCcore-8.3.0/lib;/apps/software/OpenMPI/3.1.4-gcccuda-2019b/lib;/apps/software/hwloc/1.11.12-GCCcore-8.3.0/lib;/apps/software/libpciaccess/0.14-GCCcore-8.3.0/lib;/apps/software/libxml2/2.9.9-GCCcore-8.3.0/lib;/apps/software/XZ/5.2.4-GCCcore-8.3.0/lib;/apps/software/numactl/2.0.12-GCCcore-8.3.0/lib;/apps/software/zlib/1.2.11-GCCcore-8.3.0/lib;/apps/software/CUDA/10.1.243-GCC-8.3.0/lib64/stubs;/apps/software/CUDA/10.1.243-GCC-8.3.0/lib64;/apps/software/imkl/2019.5.281-iimpi-2019b/mkl/lib/intel64;/apps/software/imkl/2019.5.281-iimpi-2019b/lib/intel64;/apps/software/impi/2018.5.288-iccifort-2019.5.281/lib64;/apps/software/iccifort/2019.5.281/compilers_and_libraries_2019.5.281/linux/tbb/lib/intel64/gcc4.4;/apps/software/iccifort/2019.5.281/compilers_and_libraries_2019.5.281/linux/mpi/intel64;/apps/software/iccifort/2019.5.281/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64;/apps/software/iccifort/2019.5.281/compilers_and_libraries_2019.5.281/linux/ipp/lib/intel64;/apps/software/iccifort/2019.5.281/compilers_and_libraries_2019.5.281/linux/daal/lib/intel64_lin;/apps/software/binutils/2.32-GCCcore-8.3.0/lib;/apps/software/GCCcore/8.3.0/lib")
set(CMAKE_C_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
