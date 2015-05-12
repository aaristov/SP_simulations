% 3D xcorr test

% 30 April 2015
% Several spot search using absolute maximum of the xcorr2 matrix
% The algorithm misses the spots if there is more than 2 on the image.

% 4 may 2015 
% introducing an individual spot search

% 5 may
% refining 3D multy search
% 

FOV=50; %FOV 50 um
zrange=2500; % range for random zp in nm
zrangeoffset=0;
% M=zeros(PSFzframes,'gpuArray'); % Vector for the maximum of xcorr2
% indM=zeros(PSFzframes,'gpuArray'); % index of maximum of xcorr2
testlength=1;
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

zstatarray=zeros(1,testlength); % (zfound)
statarray=zeros(6,testlength); % (xp,yp,zp,xerr,yerr,zerr)
zvector=-zrange/2:zrange/(testlength-1):zrange/2;
foundnum=0;
setnum=0;
maxspotnum=20;
xcorr_tr=900; % truncate mask to avoid small peaks
setarray=zeros(3,10);
foundarray=zeros(3,10);

max_C2=zeros(PSFzframes,576);
max_C2_tr=zeros(PSFzframes,576);
max_C2_tr_mask=zeros(PSFzframes,576);
imax=zeros(PSFzframes,576);
% figure(5);
% hold off;
% scatter3(1,1,1,1);
for j=1:testlength    % Number of frames with random x,y,z
    figure(5);
    hold on;
    
    [Frame,setarray]= genSPFrame(9,15,100,50,zrange,inputph,sphere,truncatecirle);
      
    figure(1);imagesc(Frame);colormap(hot);title('Original image');
%     figure(5);
    % remove noise
%     Frame=Frame-ones(size(Frame)) * mean(Frame(:));
    [thr,sorh,keepapp] = ddencmp('den','wv',Frame);
    thr 
    Frame = wdencmp('gbl',Frame,'sym4',2,thr,sorh,keepapp);
    
    Frame=Frame-ones(size(Frame)) * mean(Frame(:));
    
    figure(2);imagesc(Frame);colormap(hot);title('Denoised image');
    hold on;
    
    for i=1:PSFzframes % Guessing the appropriate PSF library frame
        PSF=(PSFarraysm(:,:,i)); % The cropped PSF frame from the library
        C2=xcorr2(Frame,PSF); % Find correlation image
        [max_C2(i,:), imax(i,:)] = max(C2);
        [ymax_C2(i,:), yimax(i,:)] = max(C2,[],2);
    end
    
    % Correction for xcorr2 slope
    max_C2=max_C2.*xcorr2normmartix2; 
    ymax_C2=ymax_C2.*xcorr2normmartix2;
    
    % Apply threshold
    threshold = .5;
    max_C2_tr_mask = max_C2 > max(max_C2(:))*threshold;
    ymax_C2_tr_mask = ymax_C2 > max(ymax_C2(:))*threshold;
    max_C2_tr = max_C2 .* max_C2_tr_mask;
    ymax_C2_tr = ymax_C2 .* ymax_C2_tr_mask;
    [PSFamp, PSFind] = max(max_C2_tr); % Serching for x-axis maxima
    [yPSFamp, yPSFind] = max(ymax_C2_tr); % Serching for y-axis maxima
    xcorr_tr=max(PSFamp)*.48;
    
    figure(3);imagesc(max_C2);
    figure(4);imagesc(ymax_C2);
    
    [max_PSFamp, ind_PSFamp] = max(PSFamp);
%     
%     figure(7),plot(PSFamp,'-');title('PSFamp (xind)');
%     figure(8),plot(yPSFamp,'-');title('yPSFamp (yind)');
    
    % Search the spots using x,z xcorr array
    while max_PSFamp > xcorr_tr 
        xind = ind_PSFamp;
        zind = PSFind(ind_PSFamp);
        yind = imax(zind,xind);
        
        if (PSFamp(ind_PSFamp+1) > 0 && PSFamp(ind_PSFamp-1) > 0)

            %         z fine search
            [xcoord,ycoord,zcorr] = finesearch(max_C2,ymax_C2,xind,yind,zind,zcortable,PSFzrange,PSFzframes);

            foundarray(:,foundnum+1)=...
                [xcoord,ycoord,...
                zcorr];
            foundnum=foundnum+1;
            figure(5);
            hold on;
            scatter3(xcoord,ycoord,zcorr,100);
            hold off;
        end
        tmpmask=ones(1,576);
        ytmpmask=ones(1,576);
%         peakwidth=9+int8(abs(zind-25)/5);
        for m=1:9
            tmpmask(1,xind+m-5)=0;
            ytmpmask(1,yind+m-5)=0;
        end
        PSFamp=tmpmask .* PSFamp;
        yPSFamp=ytmpmask .* yPSFamp;
        
%          figure(7),plot(PSFamp,'o');
%          figure(8),plot(yPSFamp,'o');
        [max_PSFamp, ind_PSFamp] = max(PSFamp);
    end
    
    [max_PSFamp, ind_PSFamp] = max(yPSFamp);
%      figure(8),plot(yPSFamp,'o');
    
    
    % Search the spots using y,z xcorr array
    while max_PSFamp > xcorr_tr
        yind = ind_PSFamp;
        zind = yPSFind(ind_PSFamp);
        xind = yimax(zind,yind);
        %         z fine search
        if (PSFamp(ind_PSFamp+1) > 0 && PSFamp(ind_PSFamp-1) > 0)
            [xcoord,ycoord,zcorr] = finesearch(max_C2,ymax_C2,xind,yind,zind,zcortable,PSFzrange,PSFzframes);

            foundarray(:,foundnum+1)=...
                [xcoord,ycoord,...
                zcorr];
            foundnum=foundnum+1;


%             figure(5);
%             hold on;
%             scatter3(horcorind2coord(xind),horcorind2coord(yind),zcorr,100);
%             hold off;
        end
        
        tmpmask=ones(1,576);
        ytmpmask=ones(1,576);
%         peakwidth=6+int8(abs(zind-25)/5);
        for m=1:9
            tmpmask(1,xind+m-5)=0;
            ytmpmask(1,yind+m-5)=0;
        end
        PSFamp=tmpmask .* PSFamp;
        yPSFamp=ytmpmask .* yPSFamp;
        
%          figure(7),plot(PSFamp,'o');
%          figure(8),plot(yPSFamp,'o');
        [max_PSFamp, ind_PSFamp] = max(yPSFamp);
    end
    hold off;
%     figure(2);imagesc(Frame);colormap(hot)

    j

end
% 
setarray
foundarray
figure(5);
scatter3(1,1,1,1);

hold on;
for i=1:foundnum
    scatter3(foundarray(1,i),foundarray(2,i),foundarray(3,i),100);
end
for i=1:length(setarray)
    scatter3(setarray(1,i),setarray(2,i),setarray(3,i),1000);
end

