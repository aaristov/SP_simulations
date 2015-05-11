% 3D xcorr test
FOV=50; %FOV 50 um
zrange=3000; % range for random zp in nm
zrangeoffset=0;
% M=zeros(PSFzframes,'gpuArray'); % Vector for the maximum of xcorr2
% indM=zeros(PSFzframes,'gpuArray'); % index of maximum of xcorr2
testlength=5;
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
foundnum=0;
setnum=0;
maxspotnum=2;
setarray=zeros(3,10);
foundarray=zeros(4,10);
figure(5);
hold off;
scatter3(1,1,1,1);
for j=1:testlength    % Number of frames with random x,y,z
    figure(5);
    hold on;
    Frame=zeros(2^n);
    for jj=1:int8(rand()*maxspotnum)+1
        
        
        zp=zrangeoffset+(rand()-1/2)*zrange; % Random zp is set in nm
        %
        xp=(rand()-1/2)*FOV*0.8; % random xp, yp is set in um.
        yp=(rand()-1/2)*FOV*0.8;

        %Call PSF drawing tool using generated coordinates
        [~,~,~,~,I2,~,~]=gaussianfft2(n,0.1,xp,yp,zp,inputph,sphere,truncatecirle);
        scatter3(xp,yp,zp,2000);
        Frame=Frame+I2;
        setnum=setnum+1;
        setarray(:,setnum)=[xp,yp,zp]
    end
    hold off;
    NoisyPhotonNum=10000;
    F=Frame/max(Frame(:));
    Frame=imnoise(F*NoisyPhotonNum*1e-12, 'poisson')*1e12;
   
    prevx=0;
    prevy=0;
    prevz=0;
    prevmax=0;
    figure(2);imagesc(Frame);colormap(hot);
    figure(5);
    
    hold on;
    trig=0;
    for i=1:PSFzframes % Guessing the appropriate PSF library frame
        PSF=(PSFarraysm(:,:,i)); % The cropped PSF frame from the library
        tmp=xcorr2(Frame,PSF); % Find correlation image
        %C2(:,:,i)=tmp;
%         imagesc(tmp);
        %[M(i),indM(i)]=max(tmp(:)); % Find maximum value of xcorr2 and index of that value indM
%         M(i)
        [max_C2(i), imax] = max(abs(tmp(:)));
        [ypeak, xpeak] = ind2sub(size(tmp),imax(1));
        if max_C2(i)> prevmax  
            trig=1;
%         else
%             trig=1;
        end
        recodredmax=1400;
        if max_C2(i)< prevmax && trig==1 && abs(xpeak-prevx)<1 && prevmax>900
            %spot found
            
            foundarray(:,foundnum+1)=...
                [horcorind2coord(prevx),horcorind2coord(prevy),...
                vertcorind2coord(prevz,PSFzframes,PSFzrange),prevmax];
            foundnum=foundnum+1;
            scatter3(horcorind2coord(prevx),horcorind2coord(prevy),vertcorind2coord(i,PSFzframes,PSFzrange),100);
            %             figure(3); plot(max_C2);figure(5)
            trig=0;
            %             recordedmax=prevmax;
            
        end
        prevx=xpeak;
        prevy=ypeak;
        prevz=i;
        prevmax=max_C2(i);
         scatter3(horcorind2coord(prevx),horcorind2coord(prevy),vertcorind2coord(i,PSFzframes,PSFzrange),max_C2(i)*.1);
            
    end
    hold off;
    figure(2);imagesc(Frame);colormap(hot)
%     figure(3);hold on; plot(max_C2); hold off;
%     spotsfound=findcentroids(Frame,SP_phase,truncatecirle);
%     numberofspots=length(spotsfound);
%     if numberofspots~=2
%         figure;
%         imagesc(Frame);
%     end
%     for k=1:numberofspots
%         x = spotsfound(k).Centroid(1);
%         y = spotsfound(k).Centroid(2);
%         cropFrame = cropimage(Frame,x,y,64);
%         
%         % Fine search
%         for i=1:PSFzframes % Guessing the appropriate PSF library frame
%             PSF=(PSFarraysm(:,:,i)); % The cropped PSF frame from the library
%             tmp=corr2(PSF,cropFrame); % Find correlation image
%             [M(i),indM(i)]=max(tmp(:)); % Find maximum value of xcorr2 and index of that value indM
%         end
%         Mfitarray=[M(zindex-3) M(zindex-2) M(zindex-1) M(zindex) M(zindex+1) M(zindex+2) M(zindex+3)];
%         Mfitscale=[PSFzscale(zindex-3) PSFzscale(zindex-2) PSFzscale(zindex-1) PSFzscale(zindex) PSFzscale(zindex+1) PSFzscale(zindex+2) PSFzscale(zindex+3)];
%         Mfit=fit(Mfitscale.',Mfitarray','poly2');
%         zfound=getfitmax(Mfit)
%     end
%     zcorr=getcorrectedz(zcorrtable,zfound);
%     zcorr=getcorrectedz(z2corrtable,zcorr);
%     zcorr=getcorrectedz(z3corrtable,zcorr);
%     
%     % xy precise search
%     
%     
%     [xfound,yfound] = zindtoxy(C2,zindex);
%     
%     tmp=size(C2);
%     framesize=tmp(1);
% %     xyfound=indM(zindex); % use 'this index to retreive the xy index (indM(zindex))
% %     [xfound,yfound]=ind2sub([tmp tmp],xyfound); %Convert the index to the subscripts
%   
%     xum = (xfound/framesize-0.5)*FOV*(framesize/2^n);
%     yum = (yfound/framesize-0.5)*FOV*(framesize/2^n);
% %     xyfound=indM(zindex); % use 'this index to retreive the xy index (indM(zindex))
% %     [xfound,yfound]=ind2sub(size(Frame),xyfound); %Convert the index to the subscripts
%     %zfound = (Mmaxfound(j)/PSFzframes-.5)*PSFzrange; % Convert z index to the z coordinate
% %     zfound = (gather(zindex)/PSFzframes-.5)*PSFzrange; % Convert z index to the z coordinate
% %     zcorr=zfound;%getcorrectedz(zcorrtable,zfound);
%   
% %     zcorr=getcorrectedz(z4corrtable,zcorr);
%     zerr=zcorr - zp; % Retreive data from GPU, convert zindex to the z coord (160 nm step)
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
%     statarray(:,j)=[xp yp zp xerr yerr zerr];
end
% 
setarray
foundarray

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
