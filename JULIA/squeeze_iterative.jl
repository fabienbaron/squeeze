include("oiplot.jl")
include("chain_average.jl")
# Squeeze iterative procedure
dataname = "../sample_data/2004-data1.oifits";

maxit = 4
npixs = [16, 32, 64, 128];
pixellations = [0.8, 0.4, 0.2, 0.1];
elements = [200, 500, 1000, 4000];
basedir = "../fake/";
nchains = 10;


for i=1:maxit
  npix = npixs[i];
  pixellation=pixellations[i];
  element = elements[i];
  mkpath(string(basedir,npix));
  outputbase = string(basedir, npix,"/reconst")
  if(i==1)
    run(`../bin/squeeze $dataname -w $npix -s $pixellation -e $element -o $outputbase -chains $nchains`)
  else
    init = string(basedir,"start",npixs[i], ".fits")
    run(`../bin/squeeze $dataname -w $npix -s $pixellation -e $element -o $outputbase -chains $nchains -i $init`)
  end

  avg, err = chain_average(outputbase)
  if (i<maxit)
    avgx2 = gridx2(avg);
    f = FITS(string(basedir,"start",npixs[i+1], ".fits"), "w"); write(f, avgx2); close(f);
  else
    f = FITS(string(basedir,"avg",npixs[i], ".fits"), "w"); write(f, avg); close(f);
    f = FITS(string(basedir,"err",npixs[i], ".fits"), "w"); write(f, avg); close(f);
    f = FITS(string(basedir,"low",npixs[i], ".fits"), "w"); write(f, avg-err); close(f);
    f = FITS(string(basedir,"hi",npixs[i], ".fits"), "w"); write(f, avg+err); close(f);
  end
  imdisp(avg)
end

readline()
