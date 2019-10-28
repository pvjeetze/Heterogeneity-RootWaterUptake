function [RdmsIMG]=CreateRdmDomain(DomainInput)
%% to open a predefined domain
IMG=imread(DomainInput); % give name of image
 
sprintf('Be patient! Randomization of flow domain may take a few minutes')
%% to zoom on polygon region
IMG=double(rgb2gray(IMG(30:end-30,3:end-3, :)));% to sharpen the edges
IMG=abs((IMG-min(min(IMG)))./(max(max(IMG))-min(min(IMG)))); % to normalize gray value between zero and one
IMG(IMG<0.001)=0.001;  % to impose the smallest value to a value bigger than zero
 
%% to sharpen the domain and remove noises
[a,b]=hist(IMG(:),unique(IMG));
IDX=sortrows([a',b], 'descend');
IMG_tem2=0;
 
for i=1:length(IDX)
 
    if sum(IMG_tem2(:)~=0)<numel(IMG_tem2)
        IMG_tem=imopen(IMG==IDX(i,2), strel('square',3));
 
        if sum(sum(IMG_tem))>0
            IMG_tem=imdilate(IMG_tem, strel('square',5));
 
            if i>1
                IMG_tem2=IMG_tem2+IDX(i,2)*(((IMG_tem2~=0)-IMG_tem)<0);
            else
                IMG_tem2=IMG_tem2+IDX(i,2)*IMG_tem;
            end
        end
    end
end
 
IMG=IMG_tem2;
clear IMG_tem2 IMG_tem OffestEnd OffestStart IDX a b
 
%% to randomize domain compartments
NrCond=10;  % number of conductivities for model experiment
RdmsIMG= ToRdmsCompartments(IMG, NrCond);
end


