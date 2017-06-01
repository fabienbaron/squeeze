# gather common display tasks
using PyPlot,PyCall

#@pyimport mpl_toolkits.axes_grid1 as axgrid


#############################################
#
# General purpose plotting
#
#############################################

# double check by plotting uv coverage

function uvplot(u,v)
  fig = figure("UV plot",figsize=(10,10),facecolor="White")
  ax = axes()
  scatter(u, v,alpha=0.5)
  scatter(-u, -v,alpha=0.5)
  title("UV coverage")
  xlabel("U")
  ylabel("V")
  grid("on")
  PyPlot.draw();PyPlot.pause(0.5); # this is used to see plots when running code in batch mode
end


function v2plot_modelvsdata(baseline_v2,v2_data,v2_data_err, v2_model) #plots V2 data vs v2 model
  fig = figure("V2 plot - Model vs Data",figsize=(10,10),facecolor="White")
  subplot(211)
  errorbar(baseline_v2,v2_data,yerr=v2_data_err,fmt="o", color="Red")
  scatter(baseline_v2, v2_model, color="Blue")
  title("V2 - Model vs data plot")
  xlabel("Baseline")
  ylabel("V2")
  subplot(212)
  scatter(baseline_v2, (v2_model - v2_data)./v2_data_err)
  xlabel("Baseline")
  ylabel("Residuals (number of sigma)")
  PyPlot.show();PyPlot.pause(1);  # this is used to see plots when running code in batch mode
end


function v2plot(baseline_v2,v2_data,v2_data_err) # plots v2 data only
  fig = figure("V2 data",figsize=(10,10),facecolor="White")
  errorbar(baseline_v2,v2_data,yerr=v2_data_err,fmt="o", color="Red")
  title("V2 data")
  xlabel("Baseline")
  ylabel("V2")
  PyPlot.show();PyPlot.pause(0.5);  # this is used to see plots when running code in batch mode
end


function t3phiplot(baseline_t3,t3phi_data,t3phi_data_err) # plots v2 data only
  fig = figure("Closure phase data",figsize=(10,10),facecolor="White")
  errorbar(baseline_t3,t3phi_data,yerr=t3phi_data_err,fmt="o", color="Red")
  title("Closure phase data")
  xlabel("Baseline")
  ylabel("Closure phase (degrees)")
  PyPlot.show();PyPlot.pause(0.5);  # this is used to see plots when running code in batch mode
end


function imdisp(image, pixellation = -1)
  fig = figure("Image",figsize=(10,10),facecolor="White")
  # if pixellation < 0 -> no pixellation entered -> do not draw in milliarcseconds
  nx=Int64(sqrt(length(image)))
  #ax = gca()
  imshow(rotl90(reshape(image,(nx,nx))), ColorMap("gist_heat"), interpolation="none"); # uses Monnier's orientation
  #divider = axgrid.make_axes_locatable(ax)
  #cax = divider[:append_axes]("right", size="5%", pad=0.05)
  #colorbar(image, cax=cax)
  PyPlot.draw();PyPlot.pause(0.05);
end



function imdisp2d(image, pixellation = -1)
  fig = figure("Image",figsize=(10,10),facecolor="White")
  # if pixellation < 0 -> no pixellation entered -> do not draw in milliarcseconds
  #ax = gca()
  imshow(rotl90(image), ColorMap("gist_heat"), interpolation="none"); # uses Monnier's orientation
  #divider = axgrid.make_axes_locatable(ax)
  #cax = divider[:append_axes]("right", size="5%", pad=0.05)
  #colorbar(image, cax=cax)
  PyPlot.draw();PyPlot.pause(0.05);
end





############################################################
#
# Imaging on spheroids
#
############################################################


function plot2Dquad(star_geometry,i) # plots the ith quad projected onto the imaging plane
  projx = star_geometry.projx;
  projy = star_geometry.projy;
  #plots the nth quad in the 2D plane, using ABCD
  # this can be used to debug lots of stuff...
  fig = figure("Test counter",figsize=(10,10),facecolor="White");
  scatter(projx[i,:], projy[i,:]);
  annotate("A", xy=[projx[i,1];projy[i,1]], xycoords="data");
  annotate("B", xy=[projx[i,2];projy[i,2]], xycoords="data");
  annotate("C", xy=[projx[i,3];projy[i,3]], xycoords="data");
  annotate("D", xy=[projx[i,4];projy[i,4]], xycoords="data");
  PyPlot.draw()
  return 1
end



function plot3d_temperature(star_geometry) # this plots the temperature map
  quads_visible = star_geometry.quads_visible;
  corners_xyz = star_geometry.corners_xyz;
  const Poly3DCollection = PyPlot.mplot3d[:art3d][:Poly3DCollection]
  fig2 = figure("Spheroid plot",figsize=(10,10),facecolor="White");
  ax = Axes3D(fig2)
  xlabel("x");
  ylabel("y");
  zlabel("z");
  grid("on")
  ax[:set_xlim]([-2.0,2.0])
  ax[:set_ylim]([-2.0,2.0])
  ax[:set_zlim]([-2.0,2.0])
  for i=1:npix
    if (quads_visible[i] > 0)
      verts = (collect(zip(corners_xyz[i, 1,:], corners_xyz[i, 2,:], corners_xyz[i, 3,:])),);
      ax[:add_collection3d](Poly3DCollection(verts));
    end
  end
  # maybe repeat the same here with another color for hidden polygons
  PyPlot.draw()
end

function plot3d_vertices(star_geometry)
center_xyz = star_geometry.center_xyz
corners_xyz = star_geometry.corners_xyz
fig = figure("Center of healpixels",figsize=(10,10),facecolor="White");
ax = gca(projection="3d");
xlabel("x");
ylabel("y");
zlabel("z");
plot3D(center_xyz[:,1],center_xyz[:,2],center_xyz[:,3], ".", color="red");
for i=1:4
  plot3D(corners_xyz[:,1, i],corners_xyz[:,2,i],corners_xyz[:,3,i], ".", color="blue");
end
PyPlot.draw()
end

function plot2d_image(star_geometry, star_map) # this plots the temperature map onto the projected 2D image plane (= observer view)
# still missing the actual intensity (includes LD)
projx = star_geometry.projx;
projy = star_geometry.projy;
nquads_visible = star_geometry.nquads_visible;
@pyimport matplotlib as mpl
patches = mpl.pymember("patches")
fig = figure("Epoch image",figsize=(10,10),facecolor="White")
ax = fig[:add_axes]([0.05,0.05,0.85,0.85])
xlabel("x");
ylabel("y");
ax[:set_xlim]([-1,1])
ax[:set_ylim]([-1,1])
for i=1:nquads_visible
#p = patches[:Polygon](hcat(projx[i,:],projy[i,:]),closed=true,edgecolor="none",color=get_cmap("jet")(i/npix),fill="true",rasterized=false)
#p = patches[:Polygon](hcat(projx[i,:],projy[i,:]),closed=true,edgecolor="black", facecolor="none",rasterized=false)
p = patches[:Polygon](hcat(projx[i,:],projy[i,:]),closed=true,edgecolor="none",color=get_cmap("jet")(star_map/sum(star_map)),fill="true",rasterized=false)
ax[:add_patch](p);
end
ax[:plot]();
PyPlot.draw()
end
