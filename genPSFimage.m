imagesize=128;


% imnoise(double(template),'salt & pepper',0.8);
% figure(1);

% imagesc(template); colorbar;
len=1000;
xImg=cell(1,len);
xCoor=cell(1,len);
for j=1:len
    fig1=zeros(imagesize+31);
    
    template1=zeros(imagesize);
    for i=1:3
        template=rand(imagesize);
        template=max(template,0.9998);
        % figure(2);
        % imagesc(template);colorbar;
        template=template-min(template(:));
        template=template/max(template(:));
        template1=template1+template;
        % figure(3);
        % imagesc(template);colorbar;

        fig1=fig1+conv2(template,PSF.array(:,:,floor(599*rand()+1)));
    end
    xImg{1,j}=fig1(15:imagesize+14,15:imagesize+14);
    xCoor{1,j}=template1;
    figure(4);
    imagesc(xImg{1,j});colorbar;
end
figure(3);
imagesc(template1);colorbar;
figure(4);
imagesc(xImg{1,j});colorbar;
%%
for i=1:length(xCoor)
    tmp=xCoor{1,i};
    xCoor1(:,i)=tmp(:);
end
