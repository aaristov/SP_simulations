% x,y,z errors distributions:

figure(10);
plot(statarray(3,:),abs(statarray(5,:)),'ro',statarray(3,:),abs(statarray(4,:)),'go'); 
xlabel('z, nm');
ylabel('x,y error, nm');
title('x,y errors over z');
figure(11);
hold on;
plot(statarray(3,:),(statarray(6,:)),'go'); 
xlabel('z, nm');
ylabel('z error, nm');
title('z errors over z');
hold off;
figure(12);
plot(statarray(1,:),abs(statarray(6,:)),'ro',statarray(2,:),abs(statarray(6,:)),'go'); 
xlabel('x,y, nm');
ylabel('z error, nm');
title('z errors over xy');
