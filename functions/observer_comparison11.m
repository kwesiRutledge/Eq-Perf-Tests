function [ results ] = observer_comparison11( varargin )
%	observer_comparison11.m
%		The objective of this experiment is to synthesize a design that achieves
%		Equalized Recovery for a given T_m and T_a on the ACC (or other) system
%		and then test the controller using noise that follows the constraints given
%		in the problem.

	%%%%%%%%%%%%%%%%%%%
	%% Manage Inputs %%
	%%%%%%%%%%%%%%%%%%%

	if nargin == 0
		verbosity	= 0;
	elseif nargin == 1
		verbosity	= varargin{1};
	elseif nargin == 3
		verbosity	= varargin{1};
		T_missing	= varargin{2};
		T_available = varargin{3};
	elseif nargin == 4
		verbosity	= varargin{1};
		T_missing	= varargin{2};
		T_available = varargin{3};
		perf_level 	= varargin{4};
	else
		error('Unacceptable number of arguments.')
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	%Default Values for observer_comparison
	if nargin < 3
		T_missing = 3;
		T_available = 3;
	end
	
	%default Value for perf_level
	if ~exist('perf_level')
		perf_level = 1;
	end

	%Using ACC System
	load('data/system_examples/acc_p.mat');

	n = size(acc.A,1);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Synthesize Controller %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%% Missing Observations for T_missing time steps
	

	%%%%%%%%%%%%%%%%%%%%
	%% Saving Results %%
	%%%%%%%%%%%%%%%%%%%%

	results.sys = acc;
	results.experim_params.T_missing = T_missing;
	results.experim_params.T_available = T_available;
	results.experim_params.perf_level = perf_level;

end