function [ F , f , dual_vars , constraints ] = create_path_constraints( varargin )
	%Description:
	%	For a given path, use the BeliefGraph to:
	%		1. Identify which mode signals can feasibly being occurring during this path
	%		2. Create Gains for this belief path
	%		3. For each word in the final belief, create reachability constraint
	%
	%Usage:
	%	[F,f,dual_vars,constraints] = BG.create_path_constraints( path_in , P_x0, P_des , P_u , obsv_gain )
	%	[F,f,dual_vars,constraints] = BG.create_path_constraints( path_in , P_x0, P_des , P_u , obsv_gain , xhat0 )
	%
	%Inputs:
	%	path_in 	- A path of the belief graph composed of sequence of belief node indices that form the path
	%				  from the initial belief to a final 'leaf' in the graph BG.
	%	P_x0 		- The initial state set of the LCSAS.
	%	obsv_gain 	- The observer gain which uses the expected output error to correct the observer state.
	%				  We assume that the simple observer's gain is constant throughout.
	%	xhat0 		- The initial estimate of the state. Initial condition of the observer.
	%				  Plan to set this to zero most of the time.
	%
	%Assumptions:
	%	This assumes that the problem is using Disturbance feedback as outlined in the ACC 2020 draft.

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	obj = varargin{1};
	path_in = varargin{2};
	P_x0 = varargin{3};
	P_des = varargin{4};
	P_u = varargin{5};
	obsv_gain = varargin{6};

	if nargin == 7
		xhat0 = varargin{7};
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	n = size(obj.lcsas.Dyn(1).A,1);
	m = size(obj.lcsas.Dyn(1).B,2);
	p = size(obj.lcsas.Dyn(1).C,1);
	wd = size(obj.lcsas.Dyn(1).B_w,2); %Assume that the size of the disturbance doesn't change for any of the included dynamics
	vd = size(obj.lcsas.Dyn(1).C_v,2);

	cg = constr_gen();

	if ~exist('xhat0')
		xhat0 = zeros(n,1);
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Identify which mode signals are Feasible %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	last_node = obj.N(path_in(end));
	feasible_mode_L = last_node.subL;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Create Gains for this belief path %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	T_i = length(path_in);
	% Feedback Variables
	F = sdpvar(m*T_i,p*T_i,'full');
	f = sdpvar(m*T_i,1,'full');

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% For each word in the final belief, create reachability constraint %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	dual_vars = {}; constraints = [];

	for word_idx = 1:length(feasible_mode_L.words)
		%Save word
		temp_word = feasible_mode_L.words{word_idx};

		%Create Constraint
		[ Pi1 , Piu , temp_constrs ] = cg.get_robust_reachability_constraints(obj.lcsas,feasible_mode_L.words{word_idx}, ...
																				P_x0,P_des,F,f, ...
																				'P_u' , P_u, ...
																				'param_type','F1',obsv_gain,xhat0);

		dual_vars{word_idx} = {Pi1,Piu};
		constraints = constraints + temp_constrs;

	end

end