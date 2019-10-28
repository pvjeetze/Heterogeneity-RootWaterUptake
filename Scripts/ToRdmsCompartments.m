function IMG_tem=ToRdmsCompartments(IMG, NrCond)
 
%% to label each polygons
IDX=unique(IMG);
IMG_tem=0*IMG;
for i=1:length(IDX)
    
    IMG_tem=IMG_tem+ (IMG==IDX(i)).*(max(IMG_tem(:))+bwlabel(IMG==IDX(i)));
end
%
%% to give a random index to each polygon
IDX=unique(IMG_tem);
r = rand(length(IDX),1);
 
IMG_tem2=IMG_tem*0;
for i=1:length(IDX)
    IMG_tem2=IMG_tem2+r(i)*(IMG_tem==IDX(i));
end
IMG_tem=IMG_tem2;
 
 
r=linspace(min(IMG_tem(:)),max(IMG_tem(:)), NrCond);
 
IMG_tem2=IMG_tem*0;
for i=1:NrCond
    
    if i==1
        IMG_tem2=IMG_tem2+r(i)*(IMG_tem<=r(i));
    else
        IMG_tem2=IMG_tem2+r(i).*(IMG_tem<=r(i)).*(IMG_tem>r(i-1));
    end
    
end
IMG_tem=IMG_tem2;
% end
 
figure(10)
subplot(1,2,1)
imagesc(IMG)
title ('Original image')
 
subplot(1,2,2)
imagesc(IMG_tem)
title ('Randomized image')
 
 
end


