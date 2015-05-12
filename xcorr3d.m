% 3D xcorr test
FOV=50; %FOV 50 um
zrange=3000; % range for random zp in nm
zrangeoffset=0;
% M=zeros(PSFzframes,'gpuArray'); % Vector for the maximum of xcorr2
% indM=zeros(PSFzframes,'gpuArray'); % index of maximum of xcorr2
testlength=11;
% zstatarray=zeros(1,testlength);
Mmaxfound=zeros(testlength,1);
M=zeros(PSFzframes,1); % Vector for the maximum of xcorr2
M1=zeros(PSFzframes,1); % Vector for the maximum of xcorr2
% Mnormcurve=mean(Mstats')';
C2=zeros(2^n+64,2^n+64,testlength); % Vectors for the maximum of xcorr2
Mstats=zeros(PSFzframes,testlength); % Vectors for the maximum of xcorr2
M1stats=zeros(PSFzframes,testlength); % Vectors for the maximum of xcorr2
indM=zeros(PSFzframes,1); % index of maximum of xcorr2
Mvect=1:PSFzframes;
Histowidth=100; % Width of the histograms
zHisto=zeros(1,Histowidth); % z,x,y histograms arrays
xHisto=zeros(1,Histowidth);
yHisto=zeros(1,Histowidth);

zstatarray=zeros(1,testlength); % (zfound)
statarray=zeros(6,testlength); % (xp,yp,zp,xerr,yerr,zerr)
% if length(Mnormcurve)~=PSFzframes
%     Mnormcurve=zeros(PSFzframes,1);
% end
% 
% rarefy=10; %miss every n-ths PSF for the coarse search



% central peak prefit
% centalpeakprefit;
zvector=-zrange/2:zrange/(testlength-1):zrange/2;

for j=1:testlength    % Number of frames with random x,y,z
%     zp=zrangeoffset+(rand()-1/2)*zrange; % Random zp is set in nm
%     zp=double(zrangeoffset+int16((rand()-1/2)*zrange/100)*100); % Random zp is set in nm
%     zp=-1200;
      zp = -PSFzrange/2 + (PSFzrange/(testlength-1))*(j-1);
%       zp=-2500;
%       zp=0;
%     if testlength>1
%         zp=zrangeoffset+(j-1)/(testlength-1)*zrange-zrange/2; % Linear zp is set in nm
%     else
%         zp=0;
%     end
%     xp=(rand()-1/2)*FOV*0.9; % random xp, yp is set in um.
%     yp=(rand()-1/2)*FOV*0.9;
    xp=0;
    yp=0;
    %Call PSF drawing tool using generated coordinates
    [t,I1,ffimage,ph,I2,fitting,width]=gaussianfft2(n,0.1,xp,yp,zp,inputph,sphere,truncatecirle);
%     I2=gpuArray(I2); % Put into GPU the camera frame
%     F= I2/sum(I2(:));
%     NoisyPhotonNum=10000;
%     Frame=imnoise(F*NoisyPhotonNum*1e-12, 'poisson')*1e12;
%      Frame=(I2); % Set camera frame
     Frame=(I2); % Set camera frame
%     Frame=(I2)/max(I2(:)); % Normalize camera frame
%     load('zcortable5may.mat')
    % Fine search
    for i=1:PSFzframes % Guessing the appropriate PSF library frame
%         PSF=gpuArray(PSFarraysm(:,:,i)); % Put into GPU the cropped PSF frame from the library
        PSF=(PSFarraysm(:,:,i)); % The cropped PSF frame from the library
        %PSF=(PSFarraysmnorm(:,:,i)); % The cropped; and normalized PSF frame from the library
        tmp=xcorr2(Frame,PSF); % Find correlation image
%         C2=gather(C2);
        C2(:,:,i)=tmp;
        [M(i),indM(i)]=max(tmp(:)); % Find maximum value of xcorr2 and index of that value indM
%         M1(i)=M(i)-5-1/40*(i-PSFzframes/2)^2;
        [max_C2(i,:), imax(i,:)] = max(tmp);
    end
    Mstats(:,j)=M; % save M profile
    [Mmax,zindex]=max(M(:)); % Find the maximum M and it's index (basically, its zindex from the set of PSFs)
    
    [xind,yind] = ind2sub(size(tmp),indM(zindex));
    xfound = horcorind2coord(xind);
    yfound = horcorind2coord(yind);
    
%     cFrame=cropimage(Frame, int8(xfound),int8(yfound),cropsize);
%     cFrame = imcrop(I2,[2^(n-1)-cropsize/2,2^(n-1)-cropsize/2,cropsize,cropsize]);
    cFrame = imcrop(I2,[xind-32-cropsize/2,yind-32-cropsize/2,cropsize,cropsize]);
    
    for i=1:51
        aux = corrcoef(cFrame,PSFarraysm(:,:,i));r(i)=aux(1,2);
    end
    [corMax(j),zind] = max(r(:));
    zfound = zind2coord(zind,PSFzrange,PSFzframes);
    zerr=zfound - zp;
    xerr = xfound - xp;
    yerr = yfound - yp;
    
    
    
    j
%     zstatarray(:,j)=zcorr;
    statarray(:,j)=[xp yp zp xerr yerr zerr];
    indarray(:,j)=[ xind yind zind];
end
figure(1);
plot(statarray(3,:),statarray(6,:),'o');
xlabel('z coordinate, nm');
ylabel('z error, nm');

figure(2);
plot(statarray(3,:),statarray(4,:),'o');
xlabel('z coordinate, um');
ylabel('x error, um');

figure(3);
plot(statarray(3,:),statarray(5,:),'o');
xlabel('z coordinate, um');
ylabel('y error, um');

% drawhistos;
% figure(7);
% plot(PSFzscale,Mstats,'o');
% 
% title('M (correlation coeff over Z');
% % figure(8);
% plot(M1stats);
% figure(9);
% plot(mean(Mstats'));
% figure(3);
% plot(zvector,zstatarray,'o');
% title('z found vs z set');
% errorGraphsPlots;
% % figure(4);
% imagesc(PSFarraysm(:,:,zindex)); daspect([1,1,1]); title('PSF found');
% figure(5);
% imagesc(Frame); daspect([1,1,1]); title('PSF set');
% figure(6);
% imagesc(xcorr2(PSFarraysm(:,:,zindex),Frame));title('correlation found');
% figure(7);
% imagesc(xcorr2(PSFarraysm(:,:,int8((zp/PSFzrange+0.5)*PSFzframes)),Frame));title('correlation needed');
% zindexneeded=int8((zp/PSFzrange+0.5)*PSFzframes)
% zindexfound=zindex
