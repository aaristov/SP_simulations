% gaussianfft(bits, width(um), phase koeff)
clear;
load('/Users/andrey/Downloads/SP_phaseMask.mat');
RGB=imread('DH_phase.png');
DH_phase = double(rgb2gray(RGB));
n=9;
%zmax=500;
%zstep=20;
%res=zeros(zmax/zstep+1,2^n);
%FWHM=zeros(zmax/zstep+1);
%z=0:zstep:zmax;
%i=1;
figure(4);
Icam=zeros(2^n,2^n);
tx = (-2^(n-1):1:2^(n-1)-1)/2^n;
ty=tx';
fc=200;
for i=1:2^n
    for j=1:2^n
        % Cyl lens
%         Cyl(i,j)=(tx(i)^2-ty(j)^2)*fc;
        Cyl(i,j)=((tx(i)-ty(j))^2+0*(tx(i)^2+ty(j)^2))*fc;
        
        
        % Helix beam
        if ty(j) == 0
            helix(i,j)=0;
        elseif ty(j) < 0
            helix(i,j)=pi-atan(tx(i)/ty(j));
        else
            helix(i,j)=pi-atan(tx(i)/ty(j));
        end
        
       
        %Airy beam
        A = 1e-6;
        B = -1e-3;
        C = B;
        di=2^(n-1); %x-shift
        dj=2^(n-1); %y-shift
        Airy(i,j) = A * ((i + j-dj)^3 + (i - (j-dj))^3) + B * (i)^2 + C * (j-dj)^2;
    end
end

%hold on;
% for z1=0:zstep:zmax
 for xp=-20:2:+20
    zp=xp*100;
    [t,I1,ffimage,ph,I2,fit,width]=gaussianfft2(n,0.1,xp,0,zp,Airy);
    Icam=Icam+I2;
 end
    tum=t.*512*.1; % t in um
%     res(i,:)=x2;
%     FWHM(i)=width;
%     i=i+1;
    zoom=22;
    xshift=0;
 
   imagesc(tum,tum,Icam);
    colormap('hot');
    xlim([-zoom+xshift,zoom+xshift]);
    ylim([-zoom,zoom]);
title('Beam Waist');
xlabel('x,um');ylabel('y (um)');
    
%     colorscheme hot;

% end

%hold off;
% figure(1);
% plot(t.*50,x.^2,'-');
% figure(2);
% plot(t,y,t,ph,'o');
%figure(3);
% %hold off;
% figure(2);
% plot(z,FWHM); %beam width over z (1D plot)
% ylim([0,2]);
% title('Beam Waist');
% xlabel('z,um');ylabel('FWHM (um)');
% 
% figure(3);
% imagesc(tum,z,res); %Slice plot of beam profile
% xlim([-1,1]);
% title('Beam Waist');
% xlabel('x,um');ylabel('z (nm)');
