
function ToWriteVariables (Output, FunctionName,Input1, Input2)

% Save ceofficient function to file.
s_file = fileread( [mfilename('fullpath'),'.m'] );
ix = strfind( s_file, 'function' );
A=s_file(ix(end):end);
if nargin==4
    B=sprintf('function %s=%s(%s,%s)', Output,FunctionName, Input1, Input2);
elseif nargin==3
    B=sprintf('function %s=%s(%s)', Output,FunctionName, Input1);
    
end
s_fcn=[B, A(25:end)] ;

ix1 = find( s_fcn == '=', 1 );
ix2 = find( s_fcn == '(', 1 );
s_fname = strtrim(s_fcn(ix1+1:ix2-1));
fid_tmp = fopen( [s_fname,'.m'], 'w' );
fprintf( fid_tmp, '%s', s_fcn );
fclose( fid_tmp );

end

%------------------------------------------------------------------------------%
function Output = FunctionName(Input1, Input2)
% D=D_function(region,state)


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