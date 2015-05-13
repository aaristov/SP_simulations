% Convert x,y,z indices to precise coordinates using correlation matrix

function [xcoord,ycoord,zcorr] = finesearch(max_C2,ymax_C2,xind,yind,zind,PSFzrange,PSFzframes)
    % max_C2 is an array [zind,xind] with values from max(xcorr2)
    
    % searching z
    % get the max_C2 line with xind
    zline=max_C2(:,xind); % 50 frames (PSFzframes) long
    if zind==1;
        zind=2;
    end
    zfitarray=[zline(zind-1) zline(zind) zline(zind+1) ];
    zfitscale=[(zind-1) (zind) (zind+1)];
    zfit=fit(zfitscale.',zfitarray','poly2');
    zindfound=getfitmax(zfit);
    zcoord=zind2coord(zindfound,PSFzrange,PSFzframes);
    zcorr=getcorrectedz(zcortable,zcoord);
%     zcorr=zcoord;
    
    
    % searching x
    % max_C2 line with zind 
    xline = max_C2(zind,:); % 576 pix long
    xfitarray=[xline(xind-1) xline(xind) xline(xind+1) ];
    xfitscale=[(xind-1) (xind) (xind+1)];
    xfit=fit(xfitscale.',xfitarray','poly2');
    xindfound=getfitmax(xfit);
    xcoord=horcorind2coord(xindfound);
    
    
    % searching y
    % ymax_C2 line with zind 
    yline = ymax_C2(zind,:); % 576 pix long
    yfitarray=[yline(yind-1) yline(yind) yline(yind+1) ];
    yfitscale=[(yind-1) (yind) (yind+1)];
    yfit=fit(yfitscale.',yfitarray','poly2');
    yindfound=getfitmax(yfit);
    ycoord=horcorind2coord(yindfound);
end