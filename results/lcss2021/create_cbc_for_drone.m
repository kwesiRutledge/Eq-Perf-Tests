% create_cbc_for_drone.m
% Description:
%	

%% Add Libraries %%

addpath(genpath('../../functions/'))

%% Constants %%

TimeHorizon = 6;
[ lcsas0 , x0 , TimeHorizon , P_target ] = get_differently_loaded_drone_lcsas('TimeHorizon',TimeHorizon);

%% Synthesis %%

[ drone_controller , info ] = lcsas0.FindConsistentBeliefController( P_target , ...
																	'SearchStrategy' , 'AscendingCardinality', ...
																	'DoOptimizationPruningWhere' , 'DuringSearch' )
save(['data/drone_data_' datestr(now,'ddmmmyyyy-HHMM') '.mat' ])