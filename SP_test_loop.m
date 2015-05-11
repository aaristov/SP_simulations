% Stat error script
testlength=10;
%figure(10);
DiffArray=zeros(testlength,Noisyzframes);
for indx=1:testlength
    SP_noisy_generation;
    PSF_corellation_2D;
    DiffArray(indx,:)=ZposDiff';
    %bar(DiffArray(index,:),'grouped');
end
figure(11);
bar(Noisyzscale,DiffArray','grouped');
xlabel('z position, nm');
ylabel('z localization error, nm');
histo;
    