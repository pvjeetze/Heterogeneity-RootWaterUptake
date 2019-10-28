
clear all
close all

MainDirectory=pwd; 
addpath([MainDirectory,'\Scripts'])
addpath([MainDirectory,'\KhEnsemble'])
DomainInput='SmallHexW51.png'; %%  name of a predefined flow domain 
%% To generate a random flow domain
[RdmsIMG]=CreateRdmDomain(DomainInput); 
%% --- Ensemble of hydraulic conductivities ---
% Load ensemble of K(h) as generated with the R-Script 'KhEnsembleVGE.R'
% van Genuchten parameters will be generated according to
% - coefficient of variation of K(h) parameters
% - number of compartments
% Get number of different conductivities
SigmaSq="6" ; % chose between geometric coefficient of variation of 0, 1, 3, 6, 12, 24
KhEnsemble = readtable(sprintf('KhEnsembleAlphaKs_SigmaSq%s_N10.csv',SigmaSq));
%% to solve flow equation
PDE_Solver(RdmsIMG, KhEnsemble, SigmaSq,MainDirectory); 

 
