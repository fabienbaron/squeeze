function congrid(im, npix)
sz=size(im)
extrapol = zeros(npix, npix)

for i=1:sz[1]
  for j=1:sz[2]
    extrapol[2*i-1:2*i,2*j-1:2*j]= im[i,j];
  end
end
return extrapol
end


# Squeeze iterative procedure
dataname = "../sample_data/2004-data1.oifits";
sizes = [16, 32, 64];
pixellations = [1.2, 0.6, 0.3];
elements = [200, 500, 1000];
basedir = "../fake/";


for i=1:3
  size = sizes[i];
  pixellation=pixellations[i];
  element = elements[i];
  outputdir = string(basedir, sizes[i], "/")
  run(`../bin/squeeze $dataname -w $size -s $pixellation -e $element -o $outputdir`)
  chain_average(dir)

end
