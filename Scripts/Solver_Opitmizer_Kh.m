
nConduct=10; 
CV_tem={'0.25', '0.5', '1', '1.5', '2'};
for i=3:5
    
%     CV=sprintf('%g', i); %"5" ; % chose between variance of 0, 0.5, 1, 1.5, 2, 2.5, 3

    CV=CV_tem{i}
%% Van genuchten parameters
VGP_tem=[0.429 0.009 0.118 1.260 0.239 0.899]; % Potting mix (Carminati et al. 2017), measured 20 November 2018 | Rep 2
HydraulicModel=0;
IDX_CV_par=[3,5,6]; %index of van genutchen parameters we need to scale them

KhEnsemble = readtable(sprintf('KhData/KhEnsembleVGP_CV%s_N%i.txt',CV, nConduct));
RangeAlpha=table2array(KhEnsemble(:,'RangeAlpha'));
RangeKs=table2array(KhEnsemble(:,'RangeKs'));
RangeLambda=table2array(KhEnsemble(:,'RangeLambda'));

X_tem=[RangeAlpha;RangeKs;RangeLambda];

Range=0.3;
lb=X_tem-sign(X_tem).*Range.*X_tem;
ub=X_tem+sign(X_tem).*Range.*X_tem;

options = optimoptions('fmincon','MaxIterations',300, 'UseParallel', 1);
[X_estimated]=fmincon(@(X_tem)Optimizer_Kh(nConduct,IDX_CV_par, VGP_tem,X_tem...,
),X_tem, [], [], [], [], lb, ub,[], options);  VGP_Fitted=X_estimated;


Optimizer_Kh(nConduct,IDX_CV_par, VGP_tem,VGP_Fitted)

VGP_Fitted_final(:,i)=VGP_Fitted; 

end

save('OptimizedVGP.mat', 'VGP_Fitted_final')
