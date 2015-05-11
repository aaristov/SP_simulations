% z coord to PSF index conversion
function zind = zcoord2ind(zcoord,PSFzrange,PSFzframes)
%zind - PSFarray frame number
PSFarrayzstep=PSFzrange/(PSFzframes-1);
zind = (zcoord+PSFzrange/2)/PSFarrayzstep+1;
end