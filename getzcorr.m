function getzcor = getcorrectedz(zcorrtable,zfound)
    zcorr = spline(zcorrtable(:,1),zcorrtable(:,2),zfound);
    getzcor=zfound-zcorr;
end