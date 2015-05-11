function [x,y] = zindtoxy(C2,zindex)
    img = C2(:,:,zindex);
    BW = img/max(img(:)) > 0.9;
    stats = regionprops((BW),'Centroid');
    x = stats.Centroid(1);
    y = stats.Centroid(2);
end
