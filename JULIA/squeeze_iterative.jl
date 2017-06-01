include("oiplot.jl")
include("chain_average.jl")
# Squeeze iterative procedure
dataname = "../sample_data/2004-data1.oifits";
npixs = [16, 32, 64];
pixellations = [1.2, 0.6, 0.3];
elements = [200, 500, 1000];
basedir = "../fake/";
nchains = 10;

i=1

for i=1:3
  npix = npixs[i];
  pixellation=pixellations[i];
  element = elements[i];
  outputbase = string(basedir, npix,"/reconst")
  if(i==1)
    run(`../bin/squeeze $dataname -w $npix -s $pixellation -e $element -o $outputbase -chains $nchains`)
  else
    init = string(basedir,"start",npixs[i], ".fits")
    run(`../bin/squeeze $dataname -w $npix -s $pixellation -e $element -o $outputbase -chains $nchains -i $init`)
  end

  avg, err = chain_average(outputbase)
  if (i<3)
    avgx2 = gridx2(avg);
    f = FITS(string(basedir,"start",npixs[i+1], ".fits"), "w"); write(f, avgx2); close(f);
  end

end
