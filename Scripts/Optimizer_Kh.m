
function [ObjFunc, X_tem]=Optimizer_Kh(nConduct,IDX_CV_par, VGP_tem,X_tem)

HydraulicModel=0;
VGP_tem2=repmat(VGP_tem, [nConduct, 1]);
VGP_tem2(:,IDX_CV_par)=reshape(X_tem, [nConduct,length(IDX_CV_par)]);
h=-linspace(10, 1e5, 2000)';


VGP_Measured=[0.428523926, 0.009293583, 0.118361388, 1.260107319, 0.238622357, 0.8992353380];

for i=1:nConduct
    
    [~,K(:,i),~] = FlowParameters(h,VGP_tem2(i,:), HydraulicModel);
    
end
%             K_Geo_Fitted=geomean(K,2);

K_Geo_Fitted= exp(sum(log(K),2)./length(K(1,:)));


[~,K_Meas,~] = FlowParameters(h,VGP_Measured, HydraulicModel);


figure(12)
drawnow
hold on
ax1=subplot(1,1,1, 'FontSize',12);
ax1.XScale='log';
ax1.YScale='log';

loglog(-h, K, 'LineWidth', 0.5)
loglog(-h, K_Meas, 'r', 'LineWidth', 2.5)
loglog(-h, K_Geo_Fitted, 'b', 'LineWidth', 2.5)
hold off
xlabel('Soil matric potential [cm]')
ylabel('Hydrualic conductivity [cm s^-^1]')
 

ObjFunc=norm(((K_Meas-K_Geo_Fitted).^2./K_Meas.^2)+(K_Meas-K_Geo_Fitted).^2./K_Geo_Fitted.^2);

end



