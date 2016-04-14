# SQUEEZE 2.6

## 1. Introduction
SQUEEZE is an image reconstruction for optical interferometry. It is designed to image complex astrophysical sources, while (optionally) modeling them simultaneously with analytic models. SQUEEZE is based on Markov Chain Monte-Carlo (MCMC) exploration of the imaging probability space, and reconstructs images and associated error bars from standard OIFITS data. SQUEEZE leverages the Open Multi-Processing (OpenMP) application programming interface to implement simulated annealing and parallel tempering, in the hope of avoiding the local minima better than classic gradient-based image reconstruction software. Another key difference is that SQUEEZE can reconstruct images using non-convex regularizers, e.g. the l0 norm for true compressed sensing.

SQUEEZE is developed by Pr Fabien Baron of Georgia State University and distributed under an open source (GPL v3) license. If you encounter bugs or if you have specific requests for additional features, models, or other enhancements, please send an email to Fabien Baron (baron@chara.gsu.edu).

Main features of SQUEEZE:

Imaging:

*    Fully polychromatic imaging with FITS output
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

SQUEEZE requires the [git](http://git-scm.com) and [cmake](http://www.cmake.org) packages to be installed on your machine.

SQUEEZE requires your compiler to be compatible with the C11 and
[OpenMP](http://openmp.org) standards. These are supported by gcc, the
Intel Compiler, the clang/LLVM compiler, but may not be natively
available on your platform yet (e.g. Mac OSX).

#### 2.1.1 Installing gcc on OSX

SQUEEZE includes all necessary libraries but an OpenMP-capable
compiler is required if you want to use parallel tempering or parallel
annealing.

Note1: though newer clang versions available for GNU/Linux do support OpenMP, OSX's default one does not at the moment.

Note2: the "gcc binary" included in OSX is actually nothing more than a hard link to the clang compiler.

We recommend you install gcc 4.9 from homebrew. You will need to install homebrew, then type:
```
brew install gcc --without-multilib
```

Before using the cmake command above, you will have to set the default
cmake compiler to be the homebrew one.
If you use the bash shell:
```
export CC=/usr/local/bin/gcc-4.9
export CXX=/usr/local/bin/g++-4.9
```
For tsch:
```
setenv CC /usr/local/bin/gcc-4.9
setenv CXX /usr/local/bin/g++-4.9
```

Note: if you are not sure which shell you have, you may type 'echo $SHELL'.

### 2.2 Installing SQUEEZE

First download the current git version of SQUEEZE using:
```
git clone https://github.com/fabienbaron/squeeze.git
```
which will create a squeeze subdirectory and download the SQUEEZE source into it.
Then install SQUEEZE by typing:
```
cd squeeze/build
cmake ..
make
```
This will configure and build both SQUEEZE's sublibraries, CFITSIO and RngStreams, then SQUEEZE itself.

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
./bin/squeeze mydata.oifits -w 64 -s 0.2 -chan 0 1.2e-6 1.35e-6 1 1.35e-6 1.43e-6 2 1.6e-6 1.8e-6
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
