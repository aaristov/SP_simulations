% 3D xcorr test
FOV=50; %FOV 50 um
zrange=3000; % range for random zp in nm
zrangeoffset=0;
% M=zeros(PSFzframes,'gpuArray'); % Vector for the maximum of xcorr2
% indM=zeros(PSFzframes,'gpuArray'); % index of maximum of xcorr2
testlength=51;
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
statarray=zeros(2,testlength); % (xp,yp,zp,xerr,yerr,zerr)
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
    end
    Mstats(:,j)=M; % save M profile
    [Mmax,zindex]=max(M(:)); % Find the maximum M and it's index (basically, its zindex from the set of PSFs)
    zfound = zind2coord(zindex,PSFzrange,PSFzframes);
    %     Mfitarray=[M(zindex-3) M(zindex-2) M(zindex-1) M(zindex) M(zindex+1) M(zindex+2) M(zindex+3)]; 
%     Mfitscale=[PSFzscale(zindex-3) PSFzscale(zindex-2) PSFzscale(zindex-1) PSFzscale(zindex) PSFzscale(zindex+1) PSFzscale(zindex+2) PSFzscale(zindex+3)];
%     Mfit=fit(Mfitscale.',Mfitarray','poly2');
%     zfound=getfitmax(Mfit);
%     zcorr=getcorrectedz(zcorrtable,zfound);
%     zcorr=getcorrectedz(z2corrtable,zcorr);
%     zcorr=getcorrectedz(z3corrtable,zcorr);
%     
    % xy precise search
    
    
%     [xfound,yfound] = zindtoxy(C2,zindex);
    
%     tmp=size(C2);
%     framesize=tmp(1);
%     xyfound=indM(zindex); % use 'this index to retreive the xy index (indM(zindex))
%     [xfound,yfound]=ind2sub([tmp tmp],xyfound); %Convert the index to the subscripts
  
%     xum = (xfound/framesize-0.5)*FOV*(framesize/2^n);
%     yum = (yfound/framesize-0.5)*FOV*(framesize/2^n);
%     xyfound=indM(zindex); % use 'this index to retreive the xy index (indM(zindex))
%     [xfound,yfound]=ind2sub(size(Frame),xyfound); %Convert the index to the subscripts
    %zfound = (Mmaxfound(j)/PSFzframes-.5)*PSFzrange; % Convert z index to the z coordinate
%     zfound = (gather(zindex)/PSFzframes-.5)*PSFzrange; % Convert z index to the z coordinate
%     zcorr=zfound;%getcorrectedz(zcorrtable,zfound);
  
%     zcorr=getcorrectedz(z4corrtable,zcorr);
    zerr=zfound - zp; % Retreive data from GPU, convert zindex to the z coord (160 nm step)
%     xerr=xum-xp;
%     yerr=yum-yp;
%     xer10 = int8(xerr/10);
%     yer10 = int8(yerr/10);
%     zer10 = int8(zerr/10);
%     xind=xer10+Histowidth/2;
%     yind=yer10+Histowidth/2;
%     zind=zer10+Histowidth/2;
%     if xind>0  && xind<Histowidth
%         xHisto(xind)=xHisto(1,xind)+1;
%     end
%     if yind>0 && yind<Histowidth
%         yHisto(yind)=yHisto(1,yind)+1;
%     end
%     if zind>0 && zind<Histowidth
%         zHisto(zind)=zHisto(1,zind)+1;
%     end
    j
%     zstatarray(:,j)=zcorr;
    statarray(:,j)=[zp zerr];
end
figure;
plot(statarray(1,:),statarray(2,:),'o');
xlabel('z coordinate, nm');
ylabel('z error, nm');

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
