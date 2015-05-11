% z PSF index z coordinate conversion
function zcoord = zind2coord(zind,PSFzrange,PSFzframes)
%zind - PSFarray frame number
PSFarrayzstep=PSFzrange/(PSFzframes-1);
zcoord = (zind-1)*PSFarrayzstep - PSFzrange/2;
end