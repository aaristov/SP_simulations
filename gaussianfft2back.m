function I2 = gaussianfft2back(I1,SLM,truncatecirle)
%gaussianfft Summary of this function goes here
%I use the gaussian source with n bits and w width
%after fourier conversion I add a quadratic phase to emulate out of focus
%in microsope and make back fourier transform to see the shape broadening.

FFE0=fft2(I1);
FPlane=abs(fftshift(FFE0));
FPlane=FPlane.*truncatecirle;
ph=fftshift(SLM);
%Back Fourier including nonlinear phase
I2=abs(ifft2(fftshift(FPlane).*exp(1i*(angle(FFE0)+ph))));

end

