% draw histos

[histox,histoy,histoz] = gethisto(statarray, 10);

figure(1);
bar(histox(:,1),histox(:,2));
title('Sigma x(x)');
xlabel('x (nm)');
ylabel('counts');

figure(2);
bar(histoy(:,1),histoy(:,2));
title('Sigma y(y)');
xlabel('y (nm)');
ylabel('counts');

figure(3);
bar(histoz(:,1),histoz(:,2));
title('Sigma z(z)');
xlabel('z (nm)');
ylabel('counts');

