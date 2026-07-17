# SQUEEZE 3.0

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

*    Parallel simulated annealing with Metropolis-Hastings moves (default)
*    Parallel tempering with Metropolis-Hastings moves (`-tempering`)

See section 3.2 for how to choose between them, and for how the annealing schedule
works.

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
This will configure and build the bundled CFITSIO, then SQUEEZE itself. The binary is
written to `./bin/squeeze`.

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
SQUEEZE reports which chain fitted best and writes it as the final image (see section 3.2.3).
*    Parallel tempering with 100 chains and full MCMC chain output
```
./bin/squeeze ./sample_data/2004-data1.oifits -w 64 -s 0.2 -chains 100 -tempering -fullchain
```
*    Polychromatic imaging, e.g. 3 channels (1.2 to 1.35 microns, 1.35-1.43 microns, and 1.6-1.8 microns).
```
./bin/squeeze mydata.oifits -w 64 -s 0.2 -wavchan 0 1.2e-6 1.35e-6 1 1.35e-6 1.43e-6 2 1.6e-6 1.8e-6
```
*    SPARCO imaging (= spectrally grey image and polychromatic modeling)

SPARCO is the default parametric model. Its six parameters are, in order:
`flux_star`, `UD` (stellar uniform-disc diameter in mas, 0 = point source),
`env_indx` (image flux power-law index), `lambda_0` (reference wavelength in m),
`flux_bg` and `bg_indx` (an over-resolved background: it dilutes the flux but
contributes no visibility). `-S` gives the corresponding step sizes; a step size of
0 holds that parameter fixed.

Here the star carries 48% of the flux at 1.6 microns and is unresolved, the image
follows lambda^-0.1, and there is no background:
```
./bin/squeeze mydata.oifits -w 64 -s 0.125 -tv 1e2 -P 0.48 0.0 -0.1 1.6e-6 0.0 -4.0 -S 0 0 0 0 0 0
```
To fit the stellar flux fraction and the environment index rather than fix them,
give those two a non-zero step size:
```
./bin/squeeze mydata.oifits -w 64 -s 0.125 -tv 1e2 -P 0.48 0.0 -0.1 1.6e-6 0.0 -4.0 -S 0.01 0 0.05 0 0 0
```

### 3.1.1 Choosing the parametric model

The model is selected at runtime with `-modeltype`; no recompilation is needed.
`squeeze -h` lists the available models with their parameters in order.

| `-modeltype` | parameters |
| --- | --- |
| `sparco` (default) | `flux_star UD env_indx lambda_0 flux_bg bg_indx` |
| `binary` | `delta_ra delta_dec flux_pri flux_sec UD_pri UD_sec frac_bw env_indx lambda_0` |

`-P` must supply exactly as many values as the chosen model takes; SQUEEZE stops with
the expected list if it does not. To add a model, follow the instructions at the top of
`src/models/modelcode.c`.

### 3.2 The minimization engines

### 3.2.1 Simulated annealing (default)

The annealing engine has **no fixed cooling schedule**. The temperature is slaved to the
current misfit: it relaxes towards

```
T  ->  chi2r / chi2_t          (chi2_t is set by -ct, default 4)
```

and is then floored at `tmin` (`-tm`, default 1). Concretely, T starts at 1, *heats up* as
it discovers how bad a random image is, and then cools only as fast as the fit improves.
There is no schedule to tune and no way to cool "too fast": the chain sets its own pace.

The consequence worth understanding is where it stops. Since `tmin` defaults to 1, once
`chi2r < chi2_t` the temperature floors at **T = 1, which is the true posterior**. The
engine then stops behaving like an optimizer and simply samples; the mean image it reports
is a posterior mean, not a best fit. That is usually what you want for imaging. If you want
posterior samples *with* error bars and evidence, parallel tempering (`-tempering`) does the
same job better and additionally gives you `logZ`.

`squeeze` also reports `Best single-frame chi2`, the single lowest-chi2 frame seen across all
chains, which is a crude MAP estimate that costs nothing.

### 3.2.2 Modes: `-tm` and `-fc` are a pair

`tmin` (`-tm`) and `chi2_t` (`-fc`) are not independent knobs. They select between two
intended modes of the engine, and **must be changed together**. The reason is that burn-in is
declared when *either*

```
T < tmin            (the schedule wants to be colder than you allow)
   or   chi2r < chi2_t     (the fit has reached its target)
```

Since the temperature settles at `chi2r/chi2_t`, the first test fires exactly when
`chi2r < chi2_t * tmin`. So `tmin` is doing double duty: temperature floor *and* convergence
criterion.

**Bayesian mean map (the default: `-tm 1 -fc 0`).** The chain cools until `chi2r < 4`, at
which point T floors at 1 (the true posterior) and the run samples it. The image you get is
a posterior mean. Convergence comes from the `T < tmin` route.

**Pseudo-maximum-entropy map (`-tm 0 -fc 1`).** Removing the floor lets the chain keep cooling
below T = 1, while `chi2_t = 1` regulates `chi2r` towards 1 instead of driving it down: the
classic "fix chi2 to its expected value and maximise entropy" reconstruction. Convergence comes
from the `chi2r < chi2_t` route. Measured on the bundled sample data, T settles around 0.6 to 0.75
and chi2r around 1:

```
./bin/squeeze ./sample_data/2004-data1.oifits -w 32 -s 0.3 -n 150 -tm 0 -fc 1
```

**Do not lower `tmin` on its own.** It removes the `T < tmin` convergence route without
supplying the other one, so no chain ever reaches burn-in: you get
`NO CHAIN REACHED BURN IN`, the final image silently falls back to the `-depth` heuristic, and
the `Best single-frame chi2` report disappears.

| flags | final T | burn-in | result |
| --- | --- | --- | --- |
| `-tm 1 -fc 0` (default) | 1.00 | reached | chi2r 0.81, posterior mean |
| `-tm 0 -fc 1` | ~0.72 | reached (via chi2_t) | chi2r 0.84, pseudo-MEM |
| `-tm 0.01` alone | 0.29 | **never** | falls back to `-depth` |

Note also that cooling harder is not a route to a better fit: at very low temperature the chain
freezes and stops exploring, settling for whatever minimum it happens to be nearest. Sampling
at T = 1 keeps moving and generally finds better optima.

### 3.2.3 Multiple chains

`-chains N` runs N **independent restarts**, not N cooperating chains: they do not share
information, and each can end up somewhere different. Which one you look at therefore matters,
especially when starting from random images (see 3.2.4).

SQUEEZE ranks the chains by the chi2 of their mean image, writes the **best** one to the
filename you gave `-o`, and prints the ranking:

```
Output -- Chain ranking (chi2r of each chain's mean image):
Output --   chain  0 : chi2r      65.7756   (rank 3 of 4)
Output --   chain  1 : chi2r       1.3882   (rank 1 of 4)  <-- best
Output --   chain  2 : chi2r      63.0660   (rank 2 of 4)
Output --   chain  3 : chi2r     618.2817   (rank 4 of 4)
Output -- Writing chain 1 (chi2r 1.3882) as the final image: result.fits
```

Every chain is still written individually as `<output>_MEAN_chainN.fits` (likewise
`_MEDIAN_` and `_MODE_`), so you can inspect the spread yourself. A large spread means the
problem is multi-modal and the "best" chain is a choice, not a consensus, so it is worth looking at
before trusting it. If even the best chain has a poor chi2r, SQUEEZE says so explicitly.

Runs are reproducible: each chain's random number stream is a pure function of its chain
index, so chain N always follows the same trajectory regardless of how the threads are
scheduled, and re-running the same command reproduces every output file byte for byte.
(Parallel tempering is not yet reproducible in this sense: there the chains also exchange
temperatures, and which pairs get tried depends on thread arrival order.)

### 3.2.4 Getting the grid right matters more than the starting image

Chains sometimes stall: instead of converging they settle at a high chi2r (hundreds) and stay
there, because the temperature is slaved to chi2 (3.2.1) and a chain that stops improving
therefore stops cooling. When this happens it is almost always a sign that the **grid is badly
chosen**, not that the run needs more iterations. A stalled chain does not improve even with
10x the iterations.

Measured on the bundled 2004 sample data, four chains from random starts (`-i randomthr`), 150
iterations:

| grid | field of view | chains reaching chi2r ~ 1 |
| --- | --- | --- |
| 32 px @ 0.3 mas | 9.6 mas | 2 of 4 |
| 32 px @ 0.2 mas | 6.4 mas | 0 of 4, the field of view no longer contains the source |
| **64 px @ 0.2 mas** | 12.8 mas | **4 of 4** (12 of 12 over three runs) |
| 128 px @ 0.1 mas | 12.8 mas | 3 of 4 |

At an appropriate grid the engine is robust and the starting image makes little difference:
both the default (all flux in one pixel) and `-i randomthr` converge every chain. Too small a
grid is what makes chains fragile. For this dataset 64x64 is the right choice, which is what
the examples in 3.1 use.

If chains do stall, the ranking output (3.2.3) makes it obvious, and the fix is a better grid
rather than more iterations. Note that cooling harder (`-ct`) does *not* rescue a stalled
chain. It freezes it in place faster.

### 3.3 Display utilities - Visualization

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
