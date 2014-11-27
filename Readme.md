SQUEEZE 2.0 - readme by Fabien Baron

Installation on OSX
-------------------

SQUEEZE includes all necessary libraries but an OpenMP-capable
compiler is required if you want to use parallel tempering or parallel
annealing. At the present time, this rules out the clang compiler
included with OSX, which OpenMP support is flaky. 

Apple also likes to redirect calls to gcc to the clang compiler, so to
check which compiler you really have on OSX, type: gcc --version

So, you have two choices on OSX:

* you may use gcc 4.8 from macports using:
sudo port install gcc48

* you can compile without OpenMP support by uncommenting the right
  line the Makefile


Usage
-----

SQUEEZE help can be invoked by typing 'squeeze -h'.


------------------------------------------

GDL/IDL utilities (squeeze_gdl.pro, squeeze_threads, plot_res.pro)
--------------

squeeze_gdl: displays the ongoing reconstruction (chi2 and regularizers,
current image, previous final image)

squeeze_threads: displays the ongoing reconstruction when using multiple
threads, e.g. when using parallel tempering

plot_res.pro: displays the final reconstructed FITS image, as well as how well
it fits the data. To be used after reconstruction.
