# SQUEEZE 2.7

## 1. Introduction
SQUEEZE is an image reconstruction for optical interferometry. It is designed to image complex astrophysical sources, while (optionally) modeling them simultaneously with analytic models. SQUEEZE is based on Markov Chain Monte-Carlo (MCMC) exploration of the imaging probability space, and reconstructs images and associated error bars from standard OIFITS data. SQUEEZE leverages the Open Multi-Processing (OpenMP) application programming interface to implement simulated annealing and parallel tempering, in the hope of avoiding the local minima better than classic gradient-based image reconstruction software. Another key difference is that SQUEEZE can reconstruct images using non-convex regularizers, e.g. the l0 norm for true compressed sensing.

SQUEEZE is developed by Pr Fabien Baron of Georgia State University and distributed under an open source (GPL v3) license. If you encounter bugs or if you have specific requests for additional features, models, or other enhancements, please send an email to Fabien Baron (fbaron@gsu.edu).

Main features of SQUEEZE:

Imaging:

*    Polychromatic imaging over a few channels (<10) with FITS output
*    Support for numerous regularizers: L0, Total Variation, Laplacian, Maximum Entropy, Dark Energy
*    Marginal likelihood computation for model selection
*    Full output of the MCMC chain including all probabilities

Modeling:
*    Uniform and limb-darkening discs, unresolved delta function, rings
*    Polychromatic support (a.k.a. SPARCO)
*    Bandwith smearing support

Minimization Engines:

*    Parallel simulated annealing with Metropolis-Hastings moves
*    Parallel tempering with Metropolis-Hastings moves

Supported data types:

*   Optical interferometric complex visibilities, differential visibilities, V2, T3 (amplitude and phase), T4 (coming soon)

## 2. Installation

### 2.1 Software requirements

SQUEEZE is designed to be cross-platform compatible. It has been
tested on several variants of GNU/Linux and on Mac OSX. Test on
additional platforms are most welcomed.

SQUEEZE requires [git](http://git-scm.com) and [cmake](http://www.cmake.org)
(version 3.20 or later) to be installed on your machine.

SQUEEZE requires your compiler to be compatible with the C11 and
[OpenMP](http://openmp.org) standards. These are supported by gcc, the
Intel Compiler, and clang/LLVM. An OpenMP-capable compiler is required if you
want to use parallel tempering or parallel annealing; without one SQUEEZE still
builds, but runs single-threaded.

#### 2.1.1 OpenMP on macOS

Apple's bundled clang does not ship the OpenMP runtime (and the "gcc" binary on
macOS is just a hard link to clang). Install a real gcc or LLVM's OpenMP via
homebrew:
```
brew install gcc          # or: brew install libomp
```
then point cmake at it, e.g.:
```
cmake -B build -DCMAKE_C_COMPILER=$(brew --prefix)/bin/gcc-14
```

### 2.2 Installing SQUEEZE

First download the current git version of SQUEEZE using:
```
git clone https://github.com/fabienbaron/squeeze.git
```
which will create a squeeze subdirectory and download the SQUEEZE source into it. Then you will need to initialize the submodule that pulls the latest OIFITSLIB (from https://github.com/jsy1001/oifitslib).
```
git submodule update --init
```

Then build SQUEEZE by typing:
```
cd squeeze
cmake -B build
cmake --build build -j
```
This will configure and build SQUEEZE's bundled libraries, CFITSIO and
RngStreams, then SQUEEZE itself. The binary is written to `./bin/squeeze`.

#### 2.2.1 Build options

| Option | Default | Meaning |
| --- | --- | --- |
| `SQUEEZE_NATIVE_ARCH` | `ON` | Tune for the build machine's CPU (`-march=native`). Set `OFF` for binaries that must run on other hardware. |
| `CMAKE_BUILD_TYPE` | `Release` | Standard CMake build types. |

For example, to build a portable binary:
```
cmake -B build -DSQUEEZE_NATIVE_ARCH=OFF
```

#### 2.2.2 A note on `-ffast-math` / `-Ofast`

SQUEEZE is deliberately **not** compiled with `-Ofast` or `-ffast-math`. Those
enable `-ffinite-math-only`, which permits the compiler to assume NaN and Inf
never occur -- and SQUEEZE relies on `isnan()`/`isinf()` to reject flagged or
corrupt V2/T3/VIS points during OIFITS import. Under `-Ofast` those checks are
folded away, bad points enter the chi2, and the likelihood (and the annealing
temperature) can go to `inf`. Please do not add these flags back.

## 3. Usage

### 3.1 Examples

Note: SQUEEZE help can be invoked by typing 'squeeze -h'.

*    Classic imaging (spectrally grey) on a 64x64 image grid, with pixel size 0.2 milli-arcseconds
```
./bin/squeeze ./sample_data/2004-data1.oifits -w 64 -s 0.2
```
*    Parallel simulated annealing with 50 chains, starting from a random image for each chain
```
./bin/squeeze ./sample_data/2004-data1.oifits -w 64 -s 0.2 -chains 50 -i randomthr
```
*    Parallel tempering with 100 chains and full MCMC chain output
```
./bin/squeeze ./sample_data/2004-data1.oifits -w 64 -s 0.2 -chains 100 -tempering -fullchain
```
*    Polychromatic imaging, e.g. 3 channels (1.2 to 1.35 microns, 1.35-1.43 microns, and 1.6-1.8 microns).
```
./bin/squeeze mydata.oifits -w 64 -s 0.2 -wavchan 0 1.2e-6 1.35e-6 1 1.35e-6 1.43e-6 2 1.6e-6 1.8e-6
```
*    SPARCO imaging (= spectrally grey image and polychromatic modeling)
```
./bin/squeeze mydata.oifits -w 64 -s 0.2 -P 1.6-e6 0.5 0.5 -2 -S 0 0.01 0.01 0.01
```

## 3.2 Display utilities - Visualization

SQUEEZE includes several visualization tools for GDL and Python (requires Astropy). With these you can:

* Follow single-chain reconstructions as they run, seeing chi2 and regularizations evolve in real time.
* Follow multi-chain reconstruction as they run, checking for chain mixing for parallel tempering or for converge for simulated annealing.
* Analyze the full MCMC probability chain of a reconstruction.
* Plot the residuals of the reconstructions.

To display and analyze the reconstruction process, a set of utilities
has been developped in several interpreted languages (IDL/GDL, PYTHON, JULIA).

* squeeze_display: displays the ongoing reconstruction (chi2 and regularizers,
current image, previous final image). This requires using the option "-temporaryfits" in SQUEEZE to continuously write chain files.

* squeeze_display_chains: displays the ongoing reconstruction when using multiple
chains, e.g. when using parallel tempering. This requires using the option "-temporaryfits" in SQUEEZE to continuously write chain files.

* plot_res: displays the final reconstructed FITS image, as well as how well
it fits the data. To be used after reconstruction.
