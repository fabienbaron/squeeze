using FITSIO
include("fbim2.jl")
# read fits file
imagefile = FITS("../sample_data/2004true.fits");
image=read(imagefile[1]);
fbim(image);
readline();
