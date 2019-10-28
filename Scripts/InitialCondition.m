function h=InitialCondition(region)
 
global IMG PixelSize VGP  HydraulicModel h0 qRoot RegionID...,
    RootConductivity ThicknessMemberan GrayValueXylem ShiftY...,
      
 
 
if sum(region.x)>0
    
    indx_X=floor((region.x)/PixelSize);
    indx_X(indx_X<=0)=1;
    
    indx_Y=floor((region.y-ShiftY)/PixelSize);
    indx_Y(indx_Y<=0)=1;
    
    VGP_tem=zeros(length(indx_X),size(VGP,2));
    %% to scale the VGP parameters
    for i=1:length(indx_X)
        ID(i)=IMG(indx_Y((i)), indx_X((i))); %#ok<AGROW>
        VGP_tem(i,:)=VGP(find(ID(i)==RegionID), :);
    end    
    %%
    
    if nargin==2   % to estimate h
        h=state.u'; %-nthroot(abs((nthroot(abs((VGP_tem(:,1) - VGP_tem(:,2))./(state.u' - VGP_tem(:,2))),(1-1./VGP_tem(:,4)))-1)), VGP_tem(:,4))./VGP_tem(:,3);
        
    else  % needed to estimate initial Theta based on h0
        h=repmat(h0, [length(indx_X), 1]);   % to calculate Theta based on h
        h(ID==GrayValueXylem)=qRoot(1)*ThicknessMemberan/RootConductivity +h0;
    end
 
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

