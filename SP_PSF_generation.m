% gaussianfft(bits, width(um), phase koeff)
%clear;
load('SP_phaseMask.mat');
% RGB=imread('DH_phase.png');
% DH_phase = double(rgb2gray(RGB))/226*2*pi;
n=9;
Icam=zeros(2^n,2^n);
tx = (-2^(n-1):1:2^(n-1)-1)/2^n;
ty=tx';
PSFzrange=5000;
PSFzframes=51;
cropsize=64;
PSFzscale=-PSFzrange/2:PSFzrange/(PSFzframes-1):PSFzrange/2;
inputph=SP_phase;
%hold on;
% for z1=0:zstep:zmax
PSFarray=zeros(2^n, 2^n, PSFzframes);
PSFarrayfft=zeros(2^n, 2^n, PSFzframes);
PSFarraysm=zeros(cropsize+1, cropsize+1, PSFzframes);
PSFarraysmnorm=zeros(cropsize+1, cropsize+1, PSFzframes);
zp=-PSFzrange/2:PSFzrange/(PSFzframes-1):PSFzrange/2;
for ind=1:PSFzframes
    [t,I1,ffimage,ph,I2,fitting,width]=gaussianfft2(n,0.1,0,0,zp(ind),inputph,sphere,truncatecirle);
    PSFarray(:,:,ind)=I2;
    %    PSFarraysm(:,:,ind)=imcrop(I2,[230,230,64,64]);
    PSFarraysm(:,:,ind)=imcrop(I2,[2^(n-1)-cropsize/2,2^(n-1)-cropsize/2,cropsize,cropsize]);
    %    PSFarrayfft(:,:,ind)=abs(fftshift(ifft2(I2)));
    tmp=PSFarraysm(:,:,ind);
    PSFarraysmnorm(:,:,ind)=tmp/max(tmp(:));
end
