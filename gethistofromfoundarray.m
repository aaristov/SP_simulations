function [histox,histoy,histoz] = gethistofromfoundarray(foundarray, step)

    xerrarray = foundarray(4,:)-foundarray(1,:);
    yerrarray = foundarray(5,:)-foundarray(2,:);
    zerrarray = foundarray(6,:)-foundarray(3,:);
    
    xmax=max(max(xerrarray));
    xmin=min(min(xerrarray));
    
    ymax=max(max(yerrarray));
    ymin=min(min(yerrarray));
    
    zmax=max(max(zerrarray));
    zmin=min(min(zerrarray));
    
    xvector=int16((xmin*1000-step:step:xmax*1000+step)/step)*step;
    yvector=int16((ymin*1000-step:step:ymax*1000+step)/step)*step;
    zvector=int16((zmin-step:step:zmax+step)/step)*step;
    
    histox=zeros(length(xvector),2);
    histoy=zeros(length(yvector),2);
    histoz=zeros(length(zvector),2);
    
    [~,z0ind]=min(abs(zvector));
    [~,x0ind]=min(abs(xvector));
    [~,y0ind]=min(abs(yvector));
    
    
    
    histox(:,1) = xvector;
    histoy(:,1) = yvector;
    histoz(:,1) = zvector;
    
    for i=1:length(foundarray(1,:))
        
        xer10 = int16(xerrarray(i)*1000/step);
        yer10 = int16(yerrarray(i)*1000/step);
        zer10 = int16(zerrarray(i)/step);
        
        xind = x0ind+xer10;
        yind = y0ind+yer10;
        zind = z0ind+zer10;
        
%         if xind>0  && xind<length(xvector)
            histox(xind,2)=histox(xind,2) + 1;
%         end
%         if yind>0 && yind<length(yvector)
            histoy(yind,2)=histoy(yind,2) + 1;
%         end
%         if zind>0 && zind<length(zvector)
            histoz(zind,2)=histoz(zind,2) + 1;
%         end
%         figure(1);
%         bar(histox(:,1),histox(:,2));
%         title('Sigma x(x)');
%         xlabel('x (nm)');
%         ylabel('counts');
%         figure(2);
%         bar(histoy(:,1),histoy(:,2));
%         title('Sigma y(y)');
%         xlabel('y (nm)');
%         ylabel('counts');
%         figure(3);
%         bar(histoz(:,1),histoz(:,2));
%         title('Sigma z(z)');
%         xlabel('z (nm)');
%         ylabel('counts');
    end
    
        
    
    