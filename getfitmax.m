function zindex = getfitmax(Mfit)
coeffs=coeffvalues(Mfit); % Mfit(x) = p1*x^2 + p2*x + p3
p1=coeffs(1);
p2=coeffs(2);
zindex=-p2/(2*p1); % Its a derivative
end
