function  stats  = findcentroids( Frame,SP_phase,truncatecirle)
%finding the centroids of the spots
%   Precess by FFT to enlarge double-lobe PSF, apply gaussian filter and
%   locate the spots
FF=gaussianfft2back(Frame,1.5*(SP_phase),truncatecirle);
myfilter = fspecial('gaussian',[10 10], 3);
FF2 = imfilter(FF, myfilter, 'replicate');
FFBW=FF2/max(FF2(:))>0.1;
stats = regionprops((FFBW),'Centroid');
end

