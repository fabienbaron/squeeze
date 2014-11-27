# SQUEEZE 2.0 

## Installation

First download the current git version of SQUEEZE using:
```ruby
git clone https://gitorious.org/squeeze/squeeze.git
```
then install SQUEEZE by typing:
```ruby
cd squeeze/build
cmake ..
make
```

### OSX

SQUEEZE includes all necessary libraries but an OpenMP-capable
compiler is required if you want to use parallel tempering or parallel
annealing. At the present time, this rules out the clang compiler included with
OSX, whose OpenMP support is flaky. Be warned that the default "gcc
binary" included in OSX is actually nothing more than a hard link to the clang compiler...

We recommend you use gcc 4.8 from macports. You will need to install Macports, then type:
```ruby
sudo port install gcc48
```

## Usage

SQUEEZE help can be invoked by typing 'squeeze -h'.



## GDL/PYTHON/JULIA utilities (squeeze_display, squeeze_threads, plot_res)

squeeze_display: displays the ongoing reconstruction (chi2 and regularizers,
current image, previous final image)

squeeze_threads: displays the ongoing reconstruction when using multiple
threads, e.g. when using parallel tempering

plot_res: displays the final reconstructed FITS image, as well as how well
it fits the data. To be used after reconstruction.
