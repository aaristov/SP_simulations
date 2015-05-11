
% Draw error histogram
dim=size(DiffArray);
xdim=dim(1);
ydim=dim(2);
Histo=zeros(1,41);
Histox=-200:10:200;
for i=1:xdim*ydim
    err=DiffArray(i);
    er10 = int8(err/10);
    Histo(1,er10+21)=Histo(1,er10+21)+1;
end

Histofit=fit(Histox.',Histo.','gauss1')
figure;
hold on;
bar(Histox,Histo);
plot(Histofit);
xlabel('z error, nm');
ylabel('Number of counts');
hold off;
