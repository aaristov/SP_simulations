% 2D Gaussian FFT example

clear
n=9; %bits
k=1; %Phase multiplier
fc =1; %cylindrical focus
fs=1;  %spherical lens focus
x = (-2^(n-1):1:2^(n-1)-1)/2^n; % Space vector
xc = 2^(-2)*0;
yc = 2^(-3)*0;
y = x'; % Space vector
we = 2^(-10); %Width
%x = exp(-t.^2/(we^2/2)).*exp(-1i*k*t.^2);   % Mag
%nlph=k*t.^2;
E0=exp(-(y-yc).^2/(we^2/2))*(exp(-(x-xc).^2/(we^2/2))); %2-D Gaussian

% Phase profile templates
G=zeros(length(x),length(y));
G1=zeros(length(x),length(y));
helix=zeros(length(x),length(y));
sp=zeros(length(x),length(y));
for i=1:length(x)
    for j=1:length(y)
        %G(i,j)=exp(-(x(i)^2+y(j)^2)/(we^2/2)); %2-D Gaussian
        G(i,j)=x(i)^2/fc-y(j)^2/fc; %a pair of cyl lenses
        %G2(i,j)=x(i)^2/f+x(i)*y(i)/f*100; %a cyl lens
        G1(i,j)=(x(i)^2+y(j)^2)*fs; % spherical lens
        
        % Helix beam
        if y(j) == 0
            helix(i,j)=0;
        elseif y(j) < 0
            helix(i,j)=pi-atan(x(i)/y(j));
        else
            helix(i,j)=pi-atan(x(i)/y(j));
        end
        
        %Saddle-point phase profile:
            % 4 gaussian peaks approx
        we1 = 2^(-2);
        c=5; %center shift 
        xc1 = 2^(-c);
        yc1 = 2^(-c);
        xc2 = -2^(-c);
        yc2 = 2^(-c);
        xc3 = 2^(-c);
        yc3 = -2^(-c);
        xc4 = -2^(-c);
        yc4 = -2^(-c);
        %xc1 = 0;
        %yc1 = 0;
        sp(i,j) = exp((-(x(i)-xc1)^2-(y(j)-yc1)^2)/(we1)^2/2) - exp((-(x(i)-xc2)^2-(y(j)-yc2)^2)/(we1)^2/2) - exp((-(x(i)-xc3)^2-(y(j)-yc3)^2)/(we1)^2/2) + exp((-(x(i)-xc4)^2-(y(j)-yc4)^2)/(we1)^2/2);
        
        %Airy beam
        A = 1e-6;
        B = -1e-3;
        C = B;
        di=2^(n-1); %x-shift
        dj=2^(n-1); %y-shift
        AiryPix = 2^n;
        % Airy(i,j) = AiryPix^3 * A * ((x(i)+.5  + y(j))^3 + (x(i)+.5 - y(j))^3) + AiryPix^2 *  B * x(i)^2 + AiryPix^2 *  C * y(i)^2;
        Airy(i,j) = A * ((i + j-dj)^3 + (i - (j-dj))^3) + B * (i)^2 + C * (j-dj)^2;
    end
end
helix=(abs(helix)-10)*4; % optimized phase profile

%Show initial intensity profile
% frame = 50;
% figure(1);imagesc(x*frame,y*frame,E0);title('focal plane');
% xlabel('x(micron)');ylabel('y (micron)'); 
% colormap('hot');

% Define input phase profile
ph=G;
ph=50*pi*sp+G1 ;
ph = Airy-G1;
%ph=0;
%ph=5*G-10*G1;
%
ph=G1;
%Truncate phase profile by circle:
for i=1:length(x)
    for j=1:length(y)
        if (i-2^(n-1))^2+(j-2^(n-1))^2 > (2^(n-2))^2
            ph(i,j) = 0;
        end
    end
end
%Truncate by rectangle
% for i=1:length(x)
%     for j=1:length(y)
%         if abs((j-2^(n-1))) > 1.5*2^(n-3);
%             ph(i,j) = 0;
%         end
%     end
% end

% Show input phase profile
figure(2);imagesc(x,y,(ph)); colorbar
title('Input phase');
xlabel('x(micron)'); ylabel('y (micron)');
colormap('hot');

% Fourier plane
H=fft2(E0.*exp(1i*0));
FPlane=fftshift(H);
%Truncate fourier profile by circle:
% for i=1:length(x)
%     for j=1:length(y)
%         if (i-2^(n-1))^2+(j-2^(n-1))^2 > (2^(n-1))^2
%             FPlane(i,j) = 0;
%         end
%     end
% end



% %Show the intensity profile in the fourier plane
% figure(3);imagesc(x,y,abs(FPlane));
% title('Fourier plane amplitude');
% xlabel('x(micron)');ylabel('y (micron)');
% colormap('hot');

% ffprofile = abs(fftshift(H))(50:150,50:150);
% imagesc(ffprofile);

% Show the phase profile in the fourier plane
%figure(4);imshow()),[-pi pi]);title('fourier plane phase profile');

%Camera projection using phase mask
CCD=ifft2(FPlane.*exp(1i*ph));

px=16*2^n; %camera pixel size
figure(5);imagesc(x*512,y*512,abs((CCD)));
title('Camera image');
xlabel('x(px)');ylabel('y (px)');
colormap('hot');

% fftB = abs(H)*exp(i*angle(H))
% G1 = ifft2(H);
% cmin = min(min(abs(G1)));
% cmax = max(max(abs(G1)));
% figure(5); imshow(abs(G1), [cmin cmax]);
% 
% AO=zeros(500,2001);
% AI=cat(2,zeros(1001,1000-b/2),ones(1001,b+1),zeros(1001,1000-b/2));
% A=cat(1,AO,AI,AO);

%figure(2);imagesc(x,y,A);title('Aperture');xlabel('x(micron)');ylabel('y (micron)');