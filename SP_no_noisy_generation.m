% gaussianfft(bits, width(um), phase koeff)
%clear;
%load('/Users/andrey/Downloads/SP_phaseMask.mat');
% RGB=imread('DH_phase.png');
% DH_phase = double(rgb2gray(RGB))/226*2*pi;
n=9;
Icam=zeros(2^n,2^n);
tx = (-2^(n-1):1:2^(n-1)-1)/2^n;
ty=tx';
Noisyzrange=4000;
Noisyzframes=50;
inputph=SP_phase;
%hold on;
% for z1=0:zstep:zmax
NoisyArray=zeros(2^n, 2^n, Noisyzframes);
Noisyzscale=-Noisyzrange/2:Noisyzrange/(Noisyzframes-1):Noisyzrange/2;
ZinitPos=zeros(1,Noisyzframes);
%ind=1;
NoisyPhotonNum=3e3;
parfor ind=1:Noisyzframes
   zp=rand()*Noisyzrange-Noisyzrange/2;
   [t,I1,ffimage,ph,I2,fitting,width]=gaussianfft2(n,0.1,0,0,zp,inputph);
   F= I2/sum(I2(:));
   LastNoisyImage=I2;%imnoise(F*NoisyPhotonNum*1e-12, 'poisson')*1e12;
   NoisyArray(:,:,ind)=LastNoisyImage;
   ZinitPos(ind)=zp;
   %ind=ind+1;
end
    %tum=t.*2^n*.1; % t in um
%     res(i,:)=x2;
%     FWHM(i)=width;
%     i=i+1;
    %zoom=22;
    %xshift=0;
    
    
% Final camera image 
% figure(1);
% imagesc(tum,tum,noisy(:,:,i-1));
% % colormap('hot');
% % xlim([-zoom+xshift,zoom+xshift]);
% % ylim([-zoom,zoom]);
% title('Beam Waist');
% xlabel('x,um');ylabel('y (um)');
    
% Phase profile
% figure(2);imagesc(tx,ty,(ph)); colorbar;
% title('Input phase');
% xlabel('x(micron)'); ylabel('y (micron)');
% colormap('hot');
