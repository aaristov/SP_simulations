% draw histos

[histox,histoy,histoz] = gethisto(statarray, 5);

histox1=histox(:,1);
histox2=histox(:,2);

histoy1=histoy(:,1);
histoy2=histoy(:,2);

histoz1=histoz(:,1);
histoz2=histoz(:,2);

figure(1);
hold off;
subplot(3,2,4);
bar(histox(:,1),histox(:,2));
title('Sigma x(x)');
xlabel('x (nm)');
ylabel('counts');

subplot(3,2,6);
bar(histoy(:,1),histoy(:,2));
title('Sigma y(y)');
xlabel('y (nm)');
ylabel('counts');

subplot(3,2,2);
bar(histoz(:,1),histoz(:,2));
title('Sigma z(z)');
xlabel('z (nm)');
ylabel('counts');

