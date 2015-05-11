function getzcor = getcorrectedz(zcortable,zfound)
    zcorr = spline(zcortable(1,:),zcortable(2,:),zfound);
    getzcor=zfound-zcorr;
end