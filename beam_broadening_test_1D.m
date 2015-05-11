% gaussianfft(bits, width(um), phase koeff)
clear;
n=9;
zmax=500;
zstep=20;
res=zeros(zmax/zstep+1,2^n);
FWHM=zeros(zmax/zstep+1);
z=0:zstep:zmax;
i=1;
figure(4);
%hold on;
for z1=0:zstep:zmax
    [t,x,y,ph,x2,fit,width]=gaussianfft(n,.01,z1);
    tum=t.*512*.1; % t in um
    res(i,:)=x2;
    FWHM(i)=width;
    i=i+1;
    plot(fit,tum,x2,'-');

end

%hold off;
% figure(1);
% plot(t.*50,x.^2,'-');
% figure(2);
% plot(t,y,t,ph,'o');
%figure(3);
%hold off;
figure(2);
plot(z,FWHM); %beam width over z (1D plot)
ylim([0,2]);
title('Beam Waist');
xlabel('z,um');ylabel('FWHM (um)');

figure(3);
imagesc(tum,z,res); %Slice plot of beam profile
xlim([-1,1]);
title('Beam Waist');
xlabel('x,um');ylabel('z (nm)');
