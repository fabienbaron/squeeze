# squeeze.c is a unity build: it #includes extract_oifits.c, regularizations.c
# and models/modelcode.c directly, so it is the only translation unit here.
add_executable(squeeze squeeze.c)

target_link_libraries(squeeze PRIVATE oitable rngstreams m)

if(OpenMP_C_FOUND)
  target_link_libraries(squeeze PRIVATE OpenMP::OpenMP_C)
endif()

# ---------------------------------------------------------------------------
# Optimisation
# ---------------------------------------------------------------------------
# Release already implies -O3. We add fp-contract (FMA) and vectorisation, but
# deliberately NOT -Ofast/-ffast-math: that enables -ffinite-math-only, under
# which the compiler assumes NaN/Inf cannot occur and folds the isnan()/isinf()
# checks in extract_oifits.c to false -- silently letting flagged OIFITS points
# into the chi2. See docs in Readme.md.
target_compile_options(squeeze PRIVATE
  $<$<C_COMPILER_ID:GNU,Clang,AppleClang,IntelLLVM>:
    -ffp-contract=fast
    -ftree-vectorize
    -ggdb
    -fno-omit-frame-pointer>)

if(SQUEEZE_NATIVE_ARCH)
  include(CheckCCompilerFlag)
  check_c_compiler_flag("-march=native" HAVE_MARCH_NATIVE)
  if(HAVE_MARCH_NATIVE)
    target_compile_options(squeeze PRIVATE -march=native)
  else()
    # Apple Silicon / some ARM toolchains use -mcpu=native instead.
    check_c_compiler_flag("-mcpu=native" HAVE_MCPU_NATIVE)
    if(HAVE_MCPU_NATIVE)
      target_compile_options(squeeze PRIVATE -mcpu=native)
    else()
      message(STATUS "Neither -march=native nor -mcpu=native supported; using generic tuning.")
    endif()
  endif()
endif()

install(TARGETS squeeze RUNTIME DESTINATION bin)
