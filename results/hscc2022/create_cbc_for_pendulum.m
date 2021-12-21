% create_cbc_for_pendulum.m
% Description:
%	

%% Add Libraries %%

addpath(genpath('../../functions/'))

%% Constants %%

TimeHorizon = 5;
[ lcsas0 , ~ , ~ , eta_v , eta_x0 , P_target ] = get_uncertain_thick_pendulum_lcsas('TimeHorizon',TimeHorizon);

%% Synthesis %%

[ pendulum_controller , info ] = lcsas0.FindConsistentBeliefController( P_target , ...
																	'SearchStrategy' , 'AscendingCardinality', ...
																	'DoOptimizationPruningWhere' , 'Nowhere' , ...
																	'RemoveBilinearityInInputConstraints', true , ...
																	'RemoveBilinearityInReachabilityConstraints', true , ...
																	'LinearizeBilinearContainment', true )
save(['data/pendulum_data_' datestr(now,'ddmmmyyyy-HHMM') '.mat' ])