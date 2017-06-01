using FITSIO
function chain_average(dir)
#get size
file = string(dir, "output_MEAN_chain0.fits");
f = FITS(file);
im = read(f[1]);
header = read_header(f[1]);
sz = size(im);
nchains = header["NCHAINS"];

fullchain = zeros(sz[1], sz[2], sz[3], nchains);

for ichain=1:nchains
  file = string(dir, "output_MEAN_chain",ichain-1,".fits");
  f = FITS(file);
  fullchain[:,:,:,ichain] = read(f[1]);
end
avg = mean(fullchain, 4);
err = std(fullchain, 4);

if(sz[3] == 1)
  avg = reshape(avg, sz[1], sz[2]);
  err = reshape(err, sz[1], sz[2]);
end

return (avg, err)
end
