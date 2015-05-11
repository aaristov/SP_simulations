% Self-correlation test for PSF dictionary

clear tmp
clear CC
for i=1:PSFzframes
tmp=(xcorr2(PSFarraysm(:,:,i),PSFarraysm(:,:,i)));
CC(:,i)=max(tmp);
end
% figure
% plot(CC)
% plot(max(CC))
selfcorprofile=max(CC);

for i=1:PSFzframes
    Mstats1(:,i)=Mstats(:,i)./selfcorprofile';
end

for i=1:PSFzframes
    [Mmax,zindex]=max(Mstats1(:,i)); % Find the maximum M and it's index (basically, its zindex from the set of PSFs)
    zfound = zind2coord(zindex,PSFzrange,PSFzframes);
    zerr1(i)=zfound-zvector(i);
end
