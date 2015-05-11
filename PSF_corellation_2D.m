%im=LastNoisyImage;
CorImPSF=zeros(PSFzframes,Noisyzframes);
ZposDiff=zeros(Noisyzframes,1);
%figure;
%hold on;
for j=1:Noisyzframes
    for i=1:PSFzframes
%         A = NoisyArray(:,:,j);
%         B = PSFarray(:,:,i);
%         
%         CorImPSF(i,j)=corr2(A,B);
%         
        
        A = gpuArray( NoisyArray(:,:,j));
        B = gpuArray(PSFarray(:,:,i));
        tmp=corr2(A,B);
        CorImPSF(i,j)=gather(tmp);
        %         aux = PSFarray(:,:,i);
        %         aux = aux(end:-1:1,end:-1:1);
        %         aux1= fft2(NoisyArray(:,:,j));
        %         aux2 = fft2(aux);
        %         aux3 = aux1.*aux2;
        %         aux4 = ifft2(aux3);
        %         aux5 = aux4.*aux4;
        %         CorImPSF(i,j) = sum(aux5(:));
    end
%     index=max(max(CorImPSF(:,j)));
%     ZposNm(j)=-Noisyzrange/2 + index * Noisyzrange;
%     ZposDiff(j)=Zpos(j)-ZposNm(j);
    
end
% figure;
% hold on;
for j=1:Noisyzframes
    CorCurve=CorImPSF(:,j)';
    CorCurvemax=max(CorCurve);
    CorCurvemin=min(CorCurve);
    CorCurveThrshold=0.5;
    NewCorCurve=CorCurve;
    for i=1:length(CorCurve)
        if NewCorCurve(i)>CorCurvemin+CorCurveThrshold*(CorCurvemax-CorCurvemin)
            NewCorCurve(i)=CorCurve(i) - CorCurveThrshold*(CorCurvemax-CorCurvemin);
        else
            NewCorCurve(i)=0;
        end
    end
    
    fitting=fit(PSFzscale.',NewCorCurve.','gauss1');
    coeffs = coeffvalues(fitting);
%     figure; 
%     plot(fitting,PSFzscale,NewCorCurve,'-');
    ZvalueFound(j)=coeffs(2);
    ZposDiff(j)=ZvalueFound(j)-ZinitPos(j);
end
% hold off;
%figure;
%bar(ZinitPos,ZposDiff);
%figure;
%plot(CorImPSF);