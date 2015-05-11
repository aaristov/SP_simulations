% Convert x,y,z indices to precise coordinates using correlation matrix
% version without z correction table to generate one

function [zcorr] = finesearchzcal(max_C2,ymax_C2,xind,yind,zind)
    % max_C2 is an array [zind,xind] with values from max(xcorr2)
    
    % searching z
    % get the max_C2 line with xind
    zline=max_C2(:,xind); % 51 frames (PSFzframes) long
    zfitarray=[zline(zind-1) zline(zind) zline(zind+1) ];
    zfitscale=[(zind-1) (zind) (zind+1)];
    zfit=fit(zfitscale.',zfitarray','poly2');
    zindfound=getfitmax(zfit);
    zcoord=zind2coord(zindfound,5000,51);
    %zcorr=getcorrectedz(zcortable,zcoord);
    zcorr=zcoord;
    % searching x
end