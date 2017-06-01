import GR

function fbim(image)
sz=size(image);
GR.clearws();
GR.setviewport(0.1, 0.9, 0.1, 0.9);
GR.setwindow(-sz[1], sz[1], -sz[1], sz[1]);
GR.setcharup(0,1)
GR.setcharheight(0.016);
GR.settextalign(2, 0);
GR.settextfontprec(3, 0);
x=collect(linspace(-sz[1],sz[1], sz[1]));
y=collect(linspace(-sz[1],sz[1], sz[1]));
GR.setspace(minimum(image), maximum(image), 0, 0);
GR.setcolormap(132);
GR.surface(x,y,image,5);
GR.textext(0.5, 0.01, "Right Ascension (mas)");
GR.setcharup(-1,0)
GR.textext(0.025, 0.5, "Declination (mas)");
GR.axes(20, 20, -256, -256, 5, 5, -0.01);
GR.updatews();
end
