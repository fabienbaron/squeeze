using FITSIO
function chain_average(outputbase)
#get size
file = string(outputbase, "_MEAN_chain0.fits");
f = FITS(file);
sz = size(read(f[1]));
header = read_header(f[1]);
nchains = header["NCHAINS"];

fullchain = zeros(sz[1], sz[2], sz[3], nchains);

for ichain=1:nchains
  file = string(outputbase, "_MEAN_chain",ichain-1,".fits");
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

function gridx2(image)
sz=size(image)

if (ndims(image) == 2)
  extrapol = zeros(2*sz[1], 2*sz[2])
  for i=1:sz[1]
    for j=1:sz[2]
      extrapol[2*i-1:2*i,2*j-1:2*j]= image[i,j];
    end
  end
else
  extrapol = zeros(2*sz[1], 2*sz[2], sz[3])
  for i=1:sz[1]
    for j=1:sz[2]
      extrapol[2*i-1:2*i,2*j-1:2*j, :]= image[i,j, :];
    end
  end
end
return extrapol
end
