function [C, K,theta] = FlowParameters(h,VGP, HydraulicModel)

% % % % script to fit van Genuchten-Mualem equation % 
% % % % Written by Dr. Mohsen Zarebanadkouki
IDX_Sat=find(h>=0);
h=-(h); 
 %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Van Genuchten; PDI  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
if HydraulicModel==0   % Van Genuchten; noraml
 
 theta_S = VGP(:,1);
    theta_R = VGP(:,2);
    alpha   = VGP(:,3);
    n       = VGP(:,4);
    m       = 1-1./n;
    Ksat    = VGP(:,5);
    lamda    = VGP(:,6);
        % Compute the effective saturation
Se=(1./(1+(alpha.*h).^n)).^m;
    % Compute the volumetric moisture content
    theta =(theta_S - theta_R).*Se + theta_R;

    % Compute the hydraulic conductivity
    K = Ksat.*Se.^(lamda).*(1 - (1 - Se.^(1./m)).^m).^2 ;

    % Compute the specific moisture storage
C = -alpha.*n.*sign(h).*(1./n - 1).*(alpha.*abs(h)).^(n - 1).*(theta_R - theta_S).*((alpha.*abs(h)).^n + 1).^(1./n - 2);

end
C=-(C); 
    

end
