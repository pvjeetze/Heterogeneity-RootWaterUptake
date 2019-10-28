
function PDE_Solver(RdmsIMG, KhEnsemble, SigmaSq,MainDirectory)

delete('K_function.m', 'C_function.m', 'Convert2Theta_Function.m')


global RdmsIMG PixelSize VGP Rroot HydraulicModel h0 qRoot RegionID...,
    RootConductivity ThicknessMemberan GrayValueXylem GrayValueMemberane ShiftY...,
    Alpha_RootMemberan Alpha_Xylem XylemConductivity
sprintf('Be patient! Numerical simulation may take a few minutes')
%% To define the size of domain;
GrayValueMemberane=1e-10; % gray value assigned to the domain representing the memberane; do not change it
GrayValueXylem=1e-20;% gray value assigned to the domain representing the xylem % do not change it
% the height will be calculated, do not change it
Rroot=5*0.02; % root radius [cm]%%%%
Width=5.1; % thickness of soil [cm]  %%%
PixelSize=Width/size(RdmsIMG,2) ;
Height=size(RdmsIMG,1)*PixelSize;   % to keep the coordinate consistent
RatioEndo2Root=0.5;
NumPixelXylem=floor((Rroot-PixelSize)*RatioEndo2Root/PixelSize);
NumPixelMemberan=floor((Rroot-PixelSize)*(1-RatioEndo2Root)/PixelSize);

%% to add root and bulk domain
% RdmsIMG(:, end-NumPixelXylem:end)= min(RdmsIMG(:));
RdmsIMG=[GrayValueMemberane*ones(size(RdmsIMG,1),NumPixelMemberan ), RdmsIMG];
RdmsIMG=[GrayValueXylem*ones(size(RdmsIMG,1),NumPixelXylem ), RdmsIMG];

% to update the thickness of entire domain
Width=size(RdmsIMG,2)*PixelSize; % thickness of soil   %%%
ThicknessMemberan=NumPixelMemberan*PixelSize;
ThicknessXylem=NumPixelXylem*PixelSize;

figure(11)
subplot(2,2,1)
imagesc([1:size(RdmsIMG,2)]*PixelSize, [1:size(RdmsIMG,1)]*PixelSize, RdmsIMG), axis tight; axis equal; colormap(bone); colorbar
title('Flow Domain')

%% to create a PDE model and generate mesh
model = createpde(1);
ShiftY=PixelSize;

Xylem = [3,4,PixelSize,ThicknessXylem,ThicknessXylem,PixelSize,  ShiftY,ShiftY,Height+ShiftY,Height+ShiftY]';
RootMemberan = [3,4,ThicknessXylem,ThicknessMemberan+ThicknessXylem,  ThicknessMemberan+ThicknessXylem,ThicknessXylem, ShiftY,ShiftY,Height+ShiftY,Height+ShiftY]';
Soil = [3,4,ThicknessMemberan+ThicknessXylem,Width,  Width,ThicknessMemberan+ThicknessXylem,  ShiftY,ShiftY,Height+ShiftY,Height+ShiftY]';

gm = [Xylem,RootMemberan,Soil];
sf = 'Xylem+RootMemberan+Soil';
ns = char('Xylem','RootMemberan', 'Soil');
ns = ns';
g = decsg(gm,sf,ns);

figure(11)
subplot(2,2,2)
geometryFromEdges(model,g);
pdegplot(model,'EdgeLabels','on', 'FaceLabels','on' )
axis tight
title('Boundary Identifier')

figure(11)
subplot(2,2,3)
Mesh=generateMesh(model, 'GeometricOrder','quadratic','Hmax',1.8e-2,'Hmin',3.6e-4, 'Hgrad',1.3 ); %#ok<NASGU>
pdeplot(model)
axis tight
title('Generate Mesh')

%% Van genuchten parameters
VGP_tem=[0.429 0.009 0.118 1.260 0.239 0.899]; % Potting mix (Carminati et al. 2017), measured 20 November 2018 | Rep 2
HydraulicModel=0;
IDX_SigmaSq_par=[3,5]; %index of varied van genuchten parameters

% % to define storage capacity of retention curve for xylem and membrane
Alpha_Xylem=1e-16;
Alpha_RootMemberan=1e-16;

% to define hydraulic conductivity of root with
%unit of [cm/s]
XylemConductivity=1; % put a big number

%$$ data are taken from table 2 in Carminati et al. 2017
RootRsistance=1.8; %kPa.s/?g
RootConductance=1/RootRsistance; %?g s^-1 kPa^-1
RootRadi=0.026; % cm
RootLength=36.6*100; % cm
FracActiveRoot=0.1;
RootConductivity = (RootConductance/(2*pi*RootRadi*RootLength*FracActiveRoot))*(10^-6/10)*ThicknessMemberan; % #ok<NOPTS> % cm/s


RegionID=unique(RdmsIMG(RdmsIMG>GrayValueMemberane));
nConduct=length(RegionID);


RangeAlpha=table2array(KhEnsemble(:,'RangeAlpha'));
RangeKs=table2array(KhEnsemble(:,'RangeKs'));
VGP=ones(length(RegionID), length(VGP_tem));
VGP(:,1)=VGP_tem(1);
VGP(:,2)=VGP_tem(2);
if sum(IDX_SigmaSq_par==3)>0;   VGP(:,3)=RangeAlpha; else; VGP(:,3)=VGP_tem(3); end
if sum(IDX_SigmaSq_par==4)>0;   VGP(:,4)=RangeN;  else; VGP(:,4)=VGP_tem(4);end
if sum(IDX_SigmaSq_par==5)>0;   VGP(:,5)=RangeKs;  else; VGP(:,5)=VGP_tem(5);end
if sum(IDX_SigmaSq_par==6)>0;   VGP(:,6)=RangeLambda;  else;   VGP(:,6)=VGP_tem(6);end


%% to add VGP of root domain
VGP=[VGP_tem; VGP_tem; VGP];
VGP(1, 3)=Alpha_Xylem;
VGP(2, 3)=Alpha_RootMemberan ;
VGP(1, 5)=XylemConductivity;
VGP(2, 5)=RootConductivity;
RegionID=[GrayValueXylem; GrayValueMemberane; RegionID];

%% % to plot the hydraulic conductivities
VGP_Fitted=[RangeAlpha;RangeKs];
Optimizer_Kh(nConduct,IDX_SigmaSq_par, VGP_tem,VGP_Fitted)

%% Maximum xylem potential
dhMax = -10000; % Maximum matric potential during model experiment

%% Initial condition
h0=-2000;  %matric potential

%% Boundary condition information
qRoot=[-1e-7 -5e-5 0];  % root water uptake rate at the root surface[cm/s];

%% To select Boundary condition IDX and Value
EdgeTop=[7,8,9];
EdgeBot=[4,5,6];
EdgeLeft=1;
EdgeRight=2;

cd([MainDirectory,'\Scripts'])
% Save periodic boundary function to file.
ToWriteVariables('K', 'K_function', 'region','state')  % output, name of function, region, puressure
SoilHydrualicConduct = @(region,state)K_function(region,state);

ToWriteVariables('C', 'C_function', 'region','state')
SoilCapacity = @(region,state)C_function(region,state);

% to convert Theta to matric potential
ToWriteVariables('Theta', 'Convert2Theta_Function', 'region','state')
cd(MainDirectory)


%% To build up the model
model.BoundaryConditions=[];
% flux at left side and constant pressure at right side and zero flux at rest BC
applyBoundaryCondition(model,'neumann','Edge',[EdgeTop, EdgeBot] ,'g', 0, 'q',0);
applyBoundaryCondition(model,'dirichlet','Edge',[EdgeRight] ,'u', h0);

% to specify the PDE coefficient
model.InitialConditions=[];
specifyCoefficients(model,'m',0,'d',SoilCapacity,'c',SoilHydrualicConduct,'a',0,'f',0);


%% Time information
% to solve it
tFinal=18000;
dt0=300;
dtShift=180; % [s] time interval for shifting flux boundary condition

% time information after maximum xylem potential is reached
dt_hXylemMax=300; % [s] time step before and after h_XylemMax is reached

model.SolverOptions.ReportStatistics = 'on';
model.SolverOptions.AbsoluteTolerance=8e-4;
model.SolverOptions.RelativeTolerance=8e-2;
model.SolverOptions.MaxIterations=1000;
model.SolverOptions.MinStep=0;


region.x=model.Mesh.Nodes(1,:);
region.y=model.Mesh.Nodes(2,:);
X_Node=model.Mesh.Nodes(1,:);
IDX_Root=find(X_Node<(Soil(3)));


%% Options for increasing flux boundary conditions

qRoot0=qRoot(1); % initial flux boundary condition

BCIncreaseStep=-5e-7;
NStepsFluxIncrease=1;
NStepsFluxDecrease=1;
hSteadyTol=10; % [cm]
hNlinTol=1000;
tRemainPlateau=1200; % [s] tRemainPlateau >= dt0

hMaxFlag=0;
NlinFlag=0;
tlist_tem=0;
IDX_hMax=1;
IDX_Nlin=2;
h_simulated=[];
h_RootSurf=[0 0];


i=1;
while  tlist_tem(end) < (tlist_tem(IDX_hMax)+tFinal)
    
    
    i=i+1;
    
    % Conditions to initiate flux boundary conditions
    
    if i==2
        ShiftStep=NStepsFluxIncrease+1;
        Flux(i)=qRoot0;
        dt=dt0;
        IDX_hEq=i;
    elseif  tlist_tem(end) >= (tlist_tem(IDX_hEq)+tRemainPlateau) && hMaxFlag~=1
        ShiftStep=1;
    end
    
    % Conditions for stepwise increase of flux boundary condition
    
    if  ShiftStep <= NStepsFluxIncrease && hMaxFlag~=1
        ShiftStep=ShiftStep+1;
        Flux(i)=Flux(end)+(BCIncreaseStep/NStepsFluxIncrease);
        dt=dtShift/NStepsFluxIncrease ;
        IDX_hEq=i;
        
    elseif i > 2 && abs(h_RootSurf(end)-h_RootSurf(end-1)) > hNlinTol && Flux(end)==Flux(end-1) && hMaxFlag~=1 || NlinFlag==1 && hMaxFlag~=1
        Flux(i)=Flux(i-1);
        dt=90;
        NlinFlag=1;
    elseif i > 2 && ShiftStep > NStepsFluxIncrease && i > 2 && hMaxFlag~=1 &&  NlinFlag==1
        Flux(i)=Flux(i-1);
        dt=90;
    elseif i > 2 && ShiftStep > NStepsFluxIncrease && i > 2 && hMaxFlag~=1
        Flux(i)=Flux(i-1);
        dt=dt0;
    end
    
    if i > 5 && abs(h_RootSurf(end)-h_RootSurf(end-1)) > hNlinTol/2 && Flux(end)==Flux(end-1) && IDX_Nlin <= 2
        IDX_Nlin=i;
    end
    
    % Conditions after h_XylemMax is reached    
    if abs(h_RootSurf(end)) > abs(h_RootSurf(IDX_Nlin-1)+dhMax) && hMaxFlag~=1
        ShiftStep=1;
        hMaxFlag=1;
        IDX_hMax=i-1;
        BCDecreaseStep=(qRoot(end)-Flux(end));
    elseif hMaxFlag~=1
        IDX_hMax=i;
    end
    
    
    
    if hMaxFlag==1 && ShiftStep <= NStepsFluxDecrease && Flux(end) <= qRoot(end)
        ShiftStep=ShiftStep+1;
        Flux(i)=Flux(end)+(BCDecreaseStep/NStepsFluxDecrease);
        dt=dt_hXylemMax/NStepsFluxDecrease;
    elseif hMaxFlag==1 && abs(h_RootSurf(end)-h_RootSurf(end-1)) >= 10*hSteadyTol
        Flux(i)=qRoot(end);
        dt=dt_hXylemMax;
    elseif hMaxFlag==1 && abs(h_RootSurf(end)-h_RootSurf(end-1)) < 10*hSteadyTol
        Flux(i)=qRoot(end);
        dt=dt0;
    end
    
    tlist_tem=[tlist_tem, tlist_tem(end)+dt];
    applyBoundaryCondition(model,'neumann','Edge',EdgeLeft,'g', Flux(i) , 'q',0);
        
    if i==2
        tlist_tem2=[tlist_tem(i-1), tlist_tem(i)];
        h0_root=Flux(i)*ThicknessMemberan/RootConductivity +h0;
        setInitialConditions(model,h0 , 'face', 3);
        setInitialConditions(model,h0_root, 'face', [1,2] );
        
        
    elseif i>2
        tlist_tem2=[tlist_tem(i-1)-tlist_tem(i-1), tlist_tem(i)-tlist_tem(i-1)];
        setInitialConditions(model,results,2);
    end
    
    results_tem=[];
    tic
    results_tem = solvepde(model,tlist_tem2);
    toc
    results=[];
    results=results_tem;
    
    if i==2
        h_simulated=results.NodalSolution;
        tlist_actual=tlist_tem(1:2);
    else
        h_simulated=[h_simulated, results.NodalSolution(:,end)]; %#ok<AGROW>
        tlist_actual=[tlist_actual, tlist_tem(i)]; %#ok<AGROW>
    end
    
    % h in xylem
    IDX_Xylem=find(X_Node<=Xylem(4));
    h_Xylem=mean(h_simulated(IDX_Xylem,:));
    
    % h at soil Surface
    h_RootSurf=h_Xylem-(Flux(1:end)/(RootConductivity/ThicknessMemberan));
    h_RootSurf(end)
    
    sprintf('solution at %g min was successful and expected final time is %g min. hMax was reached: %s.',tlist_tem(i)/60, (tlist_tem(IDX_hMax)+tFinal)/60, string(hMaxFlag==1))
    
    
end

%% Results
clear Theta h K
for i=1:length(tlist_actual)
    h_tem.u=h_simulated(:,i)';
    
    Theta(:,i)= Convert2Theta_Function(region,h_tem)';
    K(:,i)=K_function(region,h_tem)';
    
    
    % to remove the root main for visualization
    h(:,i)=h_simulated(:,i); h(IDX_Root,i)=nan;
    K(:,i)=K(:,i); K(IDX_Root,i)=nan;
    
end

h_simulated=abs(h_simulated);

%  to extract the solution for each domain

% root membrane
IDX_RootMemberan=find(X_Node<=Soil(3)); IDX_RootMemberan=IDX_RootMemberan(X_Node(IDX_RootMemberan)>=(RootMemberan(3)));
h_Memberane=mean(h_simulated(IDX_RootMemberan,:));
% xylem
IDX_Xylem=find(X_Node<=Xylem(4));
h_Xylem=mean(h_simulated(IDX_Xylem,:));
% flux
SimFlux=(-diff(mean(Theta)))'./diff(tlist_actual);
% soil Surface
h_RootSurf=h_Xylem-(abs(Flux(1:end))/(RootConductivity/ThicknessMemberan));
%hMax
hMax=abs(h_RootSurf(IDX_Nlin-1)- dhMax);
% time at hMax
t_hMax=tlist_actual(IDX_hMax)/60;
%flux at hMax
q_hMax=Flux(IDX_hMax);
% compute area after hMax is reached
hInt300=trapz(tlist_actual(IDX_hMax:end)/60, h_RootSurf(IDX_hMax:end)-2000);
% compute area until 30 Min after hMax is reached
t_PosthMax=tlist_actual(IDX_hMax+1:end)-tlist_actual(IDX_hMax);
hInt30=trapz(tlist_actual(IDX_hMax:IDX_hMax+max(find(t_PosthMax<1980)))/60,...
    h_RootSurf(IDX_hMax:IDX_hMax+max(find(t_PosthMax<1980)))-2000);
% compute area until 60 Min after hMax is reached
hInt60=trapz(tlist_actual(IDX_hMax:IDX_hMax+max(find(t_PosthMax<3780)))/60,...
    h_RootSurf(IDX_hMax:IDX_hMax+max(find(t_PosthMax<3780)))-2000);
% compute area until 180 Min after hMax is reached
hInt180=trapz(tlist_actual(IDX_hMax:IDX_hMax+max(find(t_PosthMax<10980)))/60, ...
    h_RootSurf(IDX_hMax:IDX_hMax+max(find(t_PosthMax<10980)))-2000);

%%  to print result table
table(SigmaSq, hMax,t_hMax,q_hMax,hInt30,hInt60,hInt180,hInt300)

%% to plot the results
%% Hysteresis -------------------------------------------------------------------
format short
round(max(h_RootSurf),4)
figure(13)
% Create axes
axes1 = axes();
hold(axes1,'on');
plot(-Flux(2:length(tlist_actual)),h_RootSurf(2:end),'LineStyle','-','LineWidth',1)

figure(13)
hold on
plot(-Flux(2:length(tlist_actual)),h_RootSurf(2:end),'LineStyle','-','LineWidth',1)
ylabel('\h_{root} [cm]')
xlabel('\q_{root} [cm s^{-1}]')
legend(sprintf('\sigma^2_y = %g', SigmaSq));
set(axes1,'FontName','Carlito','FontSize',15);

%% Water potential & Flux
fig=figure(14);
set(fig,'defaultAxesColorOrder',[[0 0 0]; [0 0 0.8]]);
subplot(2,3,5);
yyaxis left
plot(tlist_actual/60, h_Xylem,'Color',[0 0.5 0])
hold on
plot(tlist_actual/60, h_RootSurf,'Color',[1 0.6 0],'LineStyle','-')
hold off
ylim([15e2 2.5e4])
xlim([0 tlist_actual(end)/60])
yyaxis right
plot(tlist_actual/60, -Flux(1:length(tlist_actual)), 'Color',[0 0 0.8]);
ylim(abs([qRoot(end) -20e-6]))
%ylim([1e-7 25e-6])
yyaxis left
ylabel('Matric potential [cm]')
yyaxis right
ylabel('Simulated root water uptake [cm s^{-1}]')
xlabel('Time [min]')
legend1 = legend('Xylem potential', 'Matric potential at the root surface','Simulated root water uptake');
set(legend1,'Location','northwest','Box','off','NumColumns',2);
title(sprintf('\sigma^2_y = %g', SigmaSq))


%% Domain - matric potential ----------------------------------------------------
figure(15)
% use if you want to combine figures into one plot
% subplot(2,3,SubfigNr)
pdeplot(model,'XYData',abs(h(:,IDX_hMax)))
axis tight;
axis equal
bar=colorbar;
bar.Limits=[abs(h0), 16e+03];
caxis([abs(h0), 16e+03])
ylabel(bar, 'Soil matric potential  [cm]')
run('ColormaphMax.m')
colormap(MyCmaphMax);
title(sprintf('\sigma^2_y = %g', SigmaSq))


%% Domain - velocity-------------------------------------------------------------
results = createPDEResults(model,-h_simulated(:,IDX_hMax+2));
[cgradx,cgrady] = evaluateCGradient(results);
cgradx(IDX_Root)=nan;
cgrady(IDX_Root)=nan;

figure(16)
% use if you want to combine figures into one plot
%subplot(2,3,SubfigNr);
pdeplot(model,'XYData', (sqrt(cgradx.^2+cgrady.^2)),'ColorMap','jet','XYStyle','flat')
axis tight;
axis equal
bar=colorbar;
set(gca,'ColorScale','log')
ylabel(bar, 'Flow velocity  [cm s^-^1]')
BarMin=5e-9;
BarMax=1e-5;
caxis([BarMin  BarMax]);
bar.Limits=[BarMin  BarMax];
run('ColormapqVelocity.m')
colormap(MyCmapqVel);
title(sprintf('\sigma^2_y = %g', SigmaSq))


%% Domain - K(h) ----------------------------------------------------------------
figure(17)
% use if you want to combine figures into one plot
%subplot(2,3,SubfigNr);
K_hMax=abs(K(:,IDX_hMax));
pdeplot(model,'XYData',K_hMax, 'colorbar', 'off')
axis tight;
axis equal
bar=colorbar;
set(gca,'ColorScale','log')
BarMin=1e-12;
BarMax=4e-6;
caxis([BarMin BarMax]);
bar.Limits=[BarMin  BarMax];
ylabel(bar, 'K(h)  [cm s^-^1]')
run('ColormapKh.m')
colormap(MyCmapKh);
title(sprintf('\sigma^2_y = %g', SigmaSq))



