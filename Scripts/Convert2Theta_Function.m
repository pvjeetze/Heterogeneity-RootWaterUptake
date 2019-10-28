function Theta=Convert2Theta_Function(region,state)

global RdmsIMG PixelSize VGP  HydraulicModel  RegionID...,
    ShiftY...,
      


if sum(region.x)>0
    
    indx_X=floor((region.x)/PixelSize);
    indx_X(indx_X<=0)=1;
    
    indx_Y=floor((region.y-ShiftY)/PixelSize);
    indx_Y(indx_Y<=0)=1;
    
    VGP_tem=zeros(length(indx_X),size(VGP,2));
    %% to scale the VGP parameters
    for i=1:length(indx_X)
        ID(i)=RdmsIMG(indx_Y((i)), indx_X((i))); %#ok<AGROW>
        VGP_tem(i,:)=VGP(find(ID(i)==RegionID), :);
    end
    
    %%
    h=state.u'; 
    [C,K,Theta] = FlowParameters(h,VGP_tem, HydraulicModel);
    
    C=C';
    K=K';
    Theta=Theta';
else
    h=NaN;
    K=NaN;
    C=NaN;
    Theta=NaN;
    
end

end