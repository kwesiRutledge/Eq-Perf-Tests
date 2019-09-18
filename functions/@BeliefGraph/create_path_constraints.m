function [ Q , r , dual_vars , constraints ] = create_path_constraints(obj, path_in , P_x0, P_des , P_u)
	%Description:
	%	For a given path, use the BeliefGraph to:
	%		1. Identify which mode signals can feasibly being occurring during this path
	%		2. Create Gains for this belief path
	%		3. For each word in the final belief, create reachability constraint



	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	n = size(obj.lcsas.Dyn(1).A,1);
	m = size(obj.lcsas.Dyn(1).B,2);
	p = size(obj.lcsas.Dyn(1).C,1);
	wd = size(obj.lcsas.Dyn(1).B_w,2); %Assume that the size of the disturbance doesn't change for any of the included dynamics
	vd = size(obj.lcsas.Dyn(1).C_v,2);

	cg = constr_gen();

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
	Q = sdpvar(m*T_i,p*T_i,'full');
	r = sdpvar(m*T_i,1,'full');

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% For each word in the final belief, create reachability constraint %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	dual_vars = {}; constraints = [];

	for word_idx = 1:length(feasible_mode_L.words)
		%Save word
		temp_word = feasible_mode_L.words{word_idx};

		%Create Constraint
		[ Pi1 , Piu , temp_constrs ] = cg.get_robust_reachability_constraints(obj.lcsas,feasible_mode_L.words{word_idx}, ...
																				P_x0,P_des,Q,r, ...
																				'P_u' , P_u);

		dual_vars{word_idx} = {Pi1,Piu};
		constraints = constraints + temp_constrs;

	end

end