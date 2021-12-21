%create_cbc_for_toy3.m
%Description:
%	

%% Add Libraries %%

addpath(genpath('../../functions/'))

%% Constants %%

[ lcsas , x0 , TimeHorizon , P_target ] = get_1d_lcsas_hscc1();

%% Synthesis %%

[ toy3_controller , info ] = lcsas.FindConsistentBeliefController( P_target );

%% Visualizing %%

results.SimulationData = [];