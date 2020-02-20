function [results] = observer_comparison59( varargin )
	%Description:
	%	Testing the all mode observer construction algorithm

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	quad_sys = get_lin_quadrotor_dyn();

	n_u = size(quad_sys.Dyn(1).B,2);
	n_x = size(quad_sys.Dyn(1).A,1);
	n_y = size(quad_sys.Dyn(1).C,1);
	n_w = size(quad_sys.Dyn(1).B_w,2);
	n_v = size(quad_sys.Dyn(1).C_v,2);

	eta_u = 0.5;
	Pu = Polyhedron('lb',-eta_u*ones(1,n_u),'ub',eta_u*ones(1,n_u));

	eta_x = 1.0;
	Px0 = Polyhedron('lb',-eta_x*ones(1,n_x),'ub',eta_x*ones(1,n_x));	

	results.params.sys = quad_sys;
	results.params.Pu = Pu;
	results.params.Px0 = Px0;

	verbosity = 1; %Verbosity of Functions. Gives debugging info

	%Data file parameters.
	save_file_name = 'data/oc59_interm_results.mat';
	load_data_flag = false;
	run_opt_flag = true;

	if load_data_flag
		load(save_file_name)
		load_data_flag = true;
	end

	cg = constr_gen();

	ops = sdpsettings('verbose',verbosity);

	%%%%%%%%%%%
	%% Tests %%
	%%%%%%%%%%%

	if ~load_data_flag
		BG = BeliefGraph(quad_sys,Pu,Px0,'verbosity',verbosity);
	end
	results.BG = BG;
	
	figure;
	BG.plot()

	if run_opt_flag
		dual_var_ind = 0;
		constraints = [];

		for blang_idx = 1:length(BG.BeliefLanguage.words)
			bword_ut = BG.BeliefLanguage.words{blang_idx};

			final_lang = BG.N(bword_ut(end)).subL;
			T_i = final_lang.find_longest_length();

			%Create Q,r Gains for Each Belief Word
			Q{blang_idx} = sdpvar(n_u*T_i,n_y*T_i,'full');
			r{blang_idx} = sdpvar(n_u*T_i,1,'full');

			%Create reachability constraints for each word that could be occurring
			%(all words compatible with this belief sequence)
			for word_idx = 1:length(final_lang.words)
				possible_word = final_lang.words{word_idx};

				dual_var_ind = dual_var_ind + 1;
				[ Pi1{dual_var_ind} , Piu{dual_var_ind} , temp_constrs ] = cg.get_robust_reachability_constraints(	quad_sys, possible_word, ...
																						Px0,Px0, ...
																						Q{blang_idx},r{blang_idx}, ...
																						'P_u' , Pu);
				constraints = constraints + temp_constrs;

			end

			disp(['Created the constraints for ' num2str(blang_idx) ' belief words.' ])

		end

		%Add Block Lower Diagonal Constraints
		l_diag_constrs = cg.get_causal_constr_on_extd_gains(quad_sys,Q);
		constraints = constraints + l_diag_constrs;

		%Insert prefix constraint
		pref_constrs = cg.create_prefix_constr_on_gains( quad_sys , Q , r );
		constraints = constraints + pref_constrs;

		results.constraints = constraints;

		%Optimize!
		obj_fcn = [];
		opt_out = optimize(	constraints, obj_fcn, ops);

		%Convert Q and r matrices to F and u0
		BL = BG.BeliefLanguage;
		if opt_out.problem ~= 0
			contr = [];
		else
			% Save Feedback Matrices
			% ++++++++++++++++++++++
			Q_set = {}; r_set = {};
			F_set = {}; u0_set = {};
			for blang_idx = 1:length(BL.words)
				bword_temp = BL.words{blang_idx};
				poss_words = BG.N(bword_temp(end)).subL;

				T_i = length(poss_words.words{1});
				%Get Parameters
				[H0,S0,Cm0,J0,f_bar,B_w_big,C_v_big] = quad_sys.get_mpc_matrices('word',poss_words.words{1});

				Q_set{blang_idx} = value(Q{blang_idx});
				r_set{blang_idx} = value(r{blang_idx});
				F_set{blang_idx} = value( (inv(value(eye(size(S0,2)) + Q{blang_idx}*Cm0*S0)) ) * Q{blang_idx});
				u0_set{blang_idx} = value( inv(value(eye(size(S0,2)) + Q{blang_idx}*Cm0*S0)) * r{blang_idx} );
				
				%Fix up F and u0 to avoid NaN
				F_set{blang_idx}( isnan(F_set{blang_idx}) ) = 0;
				% u0_set{blang_idx}( isnan(u0_set{blang_idx}) ) = 0;

			end

			%Create Function Outputs
			opt_out.Q_set = Q_set;
			opt_out.r_set = r_set;

			contr = POB_Feedback(BG,F_set,u0_set);
		end

		save(save_file_name,'BG','contr')
	end

	

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Simulate the System With this Form of Controller %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	[x, u , y , sig] = contr.simulate_1run(quad_sys,Px0);

	results.x = x;
	results.u = u;
	results.y = y;
	results.sig = sig;

	figure;
	plot3(x(1,:),x(2,:),x(3,:))

	figure;
	plot(x(:,1),x(:,2))

	T = size(x,2);
	figure;
	hold on;
	plot([1:T],x(1,:))
	plot([1:T],x(2,:))
	plot([1:T],x(3,:))

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Save The Larger Variables to a Data File %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end