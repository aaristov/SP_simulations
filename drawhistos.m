% draw histos

[histox,histoy,histoz] = gethisto(statarray, 2);

histox1=histox(:,1);
histox2=histox(:,2);

histoy1=histoy(:,1);
histoy2=histoy(:,2);

histoz1=histoz(:,1);
histoz2=histoz(:,2);

histozfit = fit(histoz(:,1),histoz(:,2),'gauss1');
histoxfit = fit(histox(:,1),histox(:,2),'gauss1');
histoyfit = fit(histoy(:,1),histoy(:,2),'gauss1');

zcoefs=coeffvalues(histozfit);
ycoefs=coeffvalues(histoyfit);
xcoefs=coeffvalues(histoxfit);

zwidth=zcoefs(3);
ywidth=ycoefs(3);
xwidth=xcoefs(3);

zFWHM=2*sqrt(log(2))*zwidth;
yFWHM=2*sqrt(log(2))*ywidth;
xFWHM=2*sqrt(log(2))*xwidth;

xtitle=sprintf('Sigma x(x), FWHM %i nm',int8(xFWHM));
ytitle=sprintf('Sigma y(z), FWHM %i nm',int8(yFWHM));
ztitle=sprintf('Sigma z(z), FWHM %i nm',int8(zFWHM));
figtitle=sprintf('%d single spots test',testlength);

figure(1);
% title(figtitle);
hold off;
subplot(3,2,4);
bar(histox(:,1),histox(:,2),'w');
hold on; plot(histoxfit,'b');hold off;
title(xtitle);
xlabel('x (nm)');
ylabel('counts');

subplot(3,2,6);
bar(histoy(:,1),histoy(:,2),'w');
hold on; plot(histoyfit,'b');hold off;
title(ytitle);
xlabel('y (nm)');
ylabel('counts');

subplot(3,2,2);
bar(histoz(:,1),histoz(:,2),'w');
hold on; plot(histozfit,'b');hold off;
title(ztitle);
xlabel('z (nm)');
ylabel('counts');

