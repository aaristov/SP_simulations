% Generate SP camera frame
% spotNum - number of spots,
% photonNum - Poisson noise parameter
% gNoise - Gaussian noise ration (0 - 100)%
% n - number of bits (2^n)

function [Frame,setarray] = genSPFrame(n,spotNum,photonNum,gNoise,zrange,inputph,sphere,truncatecirle);
    

    FOV=50; % Field of view
    Frame=zeros(2^n);
    setnum=0;
    
    % Generating the camera image
%     for jj=1:int8(rand()*maxspotnum)+1
    for jj=1:spotNum
               
        zp=(rand()-1/2)*zrange; % Random zp is set in nm
        xp=(rand()-1/2)*FOV*0.8; % random xp, yp is set in um.
        yp=(rand()-1/2)*FOV*0.8;

        %Call PSF drawing tool using generated coordinates
        [~,~,~,~,I2,~,~]=gaussianfft2(n,0.1,xp,yp,zp,inputph,sphere,truncatecirle);
%         figure(5);scatter3(xp,yp,zp,1000);
        Frame=Frame+I2;
        setnum=setnum+1;
        setarray(:,setnum)=[xp,yp,zp];
    end
    
    % Add Poisson noise
    
    NoisyPhotonNum=photonNum;
    F=Frame/max(Frame(:));
    Frame=imnoise(F*NoisyPhotonNum*1e-12, 'poisson')*1e12;
    
    %Add gaussian noise
    noisepower=gNoise; % percent from frame max amplitude
    gnoise=ones(size(Frame));
    ggnoise=imnoise(gnoise,'gauss');
    Frame = Frame + ggnoise*max(Frame(:))*noisepower/100;
    
end