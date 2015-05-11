function zcorrtable = savezspline(statarray)
    a=statarray(3,:); % z set position
    b=statarray(6,:); % z found error
    zcorrtable(:,1)=a;
    zcorrtable(:,2)=b;
end