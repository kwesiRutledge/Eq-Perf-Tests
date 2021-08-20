% create_cbc_for_drone.m
% Description:
%	

%% Add Libraries %%

addpath(genpath('../../functions/'))

%% Constants %%

TimeHorizon = 4;
[ lcsas0 , x0 , TimeHorizon , P_target ] = get_differently_loaded_drone_lcsas('TimeHorizon',TimeHorizon);

%% Synthesis %%

[ drone_controller , info ] = lcsas0.FindConsistentBeliefController( P_target , 'SearchStrategy' , 'DescendingCardinality+PreferAllModeSequence' )