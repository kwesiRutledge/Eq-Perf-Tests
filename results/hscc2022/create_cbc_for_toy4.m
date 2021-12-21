%create_cbc_for_toy4.m
%Description:
%	System hscc2 example 2 in hscc.

%% Add Libraries %%

addpath(genpath('../../functions/'))

%% Constants %%

[ lcsas , x0 , TimeHorizon , P_target ] = get_1d_lcsas_hscc2();

%% Synthesis %%

[ toy4_controller , info ] = lcsas.FindConsistentBeliefController( P_target , ...
																	'RemoveBilinearityInInputConstraints', true , ...
																	'RemoveBilinearityInReachabilityConstraints', true , ...
																	'LinearizeBilinearContainment', true );

%% Visualizing %%

results.SimulationData = [];

%% Saving Data %%

save(['data/toy4_data_' datestr(now,'ddmmmyyyy-HHMM') '.mat' ])