% central peak prefit
Mpr=zeros(PSFzframes,1); % Vector for the maximum of xcorr2
Mprstats=zeros(PSFzframes,2); % Vectors for the maximum of xcorr2
M1pr=Mpr;
M2pr=Mpr;
M1prstats=zeros(PSFzframes,2); % Vectors for the maximum of xcorr2
indMpr=zeros(PSFzframes,1); % index of maximum of xcorr2
Mprvect=1:PSFzframes;
for j=1:2    % Number of frames with random x,y,z
    %zp=zrangeoffset+rand()*zrange-zrange/2; % Random zp is set in nm
    zp=(j-1.5)*.8*PSFzrange; % Linear zp is set in nm

    %Call PSF drawing tool using generated coordinates
    [t,I1,ffimage,ph,I2,fitting,width]=gaussianfft2(n,0.1,0,0,zp,inputph,sphere,truncatecirle);
%     Frame=gpuArray(I2); % Put into GPU the camera frame
    Frame=(I2)/max(I2(:)); % Put into GPU the camera frame
    
    % Fine search
    for i=1:PSFzframes % Guessing the appropriate PSF library frame
%         PSF=gpuArray(PSFarraysm(:,:,i)); % Put into GPU the cropped PSF frame from the library
        PSF=(PSFarraysmnorm(:,:,i)); % Put into GPU the cropped PSF frame from the library
        C2=xcorr2(PSF,Frame); % Find correlation image
        %C3=gather(C2);
        [Mpr(i),indMpr(i)]=max(C2(:)); % Find maximum value of xcorr2 and index of that value indM
%         M1(i)=M(i)-5-1/40*(i-PSFzframes/2)^2;
    end
    Mprstats(:,j)=Mpr; % save M profile
    %M1=M-5-1/100*(Mvect.-PSFzframes/2)^2;
    %Mprfit=fit(Mprvect.',Mpr,'poly2');
  
    
end
mean12=mean(Mprstats');
Mprfit=fit(Mprvect.',mean12','poly2');
for i=1:PSFzframes
    M1pr(i)=mean12(i)-Mprfit(i);
end
prpeakwidth= 0.15; % percentage of the whole PSFzrange

Mprmask=zeros(PSFzframes,1);
for i=PSFzframes*(1/2-prpeakwidth/2)-.5:PSFzframes*(1/2+prpeakwidth/2)-.5
    Mprmask(int16(i))=1;
end
M2pr=M1pr.*Mprmask;
figure(9);
plot(M2pr);

