function [histox,histoy,histoz] = gethisto(statarray, step)
    xmax=max(max(statarray(4,:)));
    xmin=min(min(statarray(4,:)));
    
    ymax=max(max(statarray(5,:)));
    ymin=min(min(statarray(5,:)));
    
    zmax=max(max(statarray(6,:)));
    zmin=min(min(statarray(6,:)));
    
    xvector=int16((xmin*1000:step:xmax*1000)/step)*step;
    yvector=int16((ymin*1000:step:ymax*1000)/step)*step;
    zvector=int16((zmin:step:zmax)/step)*step;
    
    histox=zeros(length(xvector),2);
    histoy=zeros(length(yvector),2);
    histoz=zeros(length(zvector),2);
    
    [mmm,z0ind]=min(abs(zvector));
    [mmm,x0ind]=min(abs(xvector));
    [mmm,y0ind]=min(abs(yvector));
    
    for i=1:length(statarray(1,:))
        
        xer10 = int16(statarray(4,i)*1000/step);
        yer10 = int16(statarray(5,i)*1000/step);
        zer10 = int16(statarray(6,i)/step);
        
        xind = x0ind+xer10;
        yind = y0ind+yer10;
        zind = z0ind+zer10;
        
        if xind>0  && xind<length(xvector)
            histox(xind,2)=histox(xind,2) + 1;
        end
        if yind>0 && yind<length(yvector)
            histoy(yind,2)=histoy(yind,2) + 1;
        end
        if zind>0 && zind<length(zvector)
            histoz(zind,2)=histoz(zind,2) + 1;
        end
        
    end
    
    histox(:,1) = xvector;
    histoy(:,1) = yvector;
    histoz(:,1) = zvector;
    
        
    
    