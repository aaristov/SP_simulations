function [t,I0,ff,nlph,I2,gss,width ] = gaussianfft( n,d,z )
%gaussianfft Summary of this function goes here
%I use the gaussian source with n bits and w width
%after fourier conversion I add a quadratic phase to emulate out of focus
%in microsope and make back fourier transform to see the shape broadening.
% The output will be the fwhm of the resulting gaussian profile.
%   Detailed explanation goes here
% n - number of bits
% d - source diameter in um
% z - z-distance out of focus
FOV=50; % Field of view of the objective (um)
t = (-2^(n-1):1:2^(n-1)-1)/2^n; % Time vector
x = exp(-t.^2/((d/FOV)^2));   % Mag
I0=x.^2; %intensity
k = z/5; %derived from Rayleigh length for wl=600nm, w0=170nm.
nlph=k*t.^2; % Non-linear phase
y=fft(x);
ff=abs(fftshift(y));
% truncate
for i=1:2^n
    if abs(i-2^(n-1)-1)>2^(n-2)
        ff(i)=0;
        nlph(i)=0;
    end
end

ph=fftshift(nlph);
%Back Fourier including nonlinear phase
x2=ifft(fftshift(ff).*exp(1i*(angle(y)+ph))); 
amp=abs(x2);
I2=amp.^2;
tum=t.*512*.1; % t in um
gss=fit(tum.',I2.','gauss1');
gssc=coeffvalues(gss);
width=gssc(3)*2*sqrt(2*log(2)); %width in nm
end

