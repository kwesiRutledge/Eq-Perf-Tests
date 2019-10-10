function [BG,contr,opt_out] = synth_robust_reach_contr( varargin )
	%Description:
	%	Testing the all mode observer construction algorithm
	%
	%Usage:
	%	[BG,contr] = lcsas.synth_robust_reach_contr( P_x0 , P_u , P_target )
	%	[BG,contr] = lcsas.synth_robust_reach_contr( P_x0 , P_u , P_target , 'BG' , BG )
	%	[BG,contr] = lcsas.synth_robust_reach_contr( P_x0 , P_u , P_target , 'debug' , debug_flag )

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	in_lcsas = varargin{1};
	P_x0 = varargin{2};
	P_u = varargin{3};
	P_target = varargin{4};

	if ~isa(in_lcsas,'LCSAS')
		error('Expected first input to be an LCSAS object.')
	end

	if ~(isa(P_x0,'Polyhedron') && isa(P_u,'Polyhedron') && isa(P_target,'Polyhedron'))
		error('Expected all sets given as input to be Polyhedron objects.')
	end

	varargin_idx = 5;
	while varargin_idx <= nargin
		switch varargin{varargin_idx}
			case 'BG'
				BG = varargin{varargin_idx+1};
				varargin_idx = varargin_idx + 2;
			case 'debug'
				debug_flag = varargin{varargin_idx+1};
				varargin_idx = varargin_idx + 2;
			otherwise
				error(['Unexpected flag for this function: ' varargin{varargin_idx}])
		end
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	n_u = size(in_lcsas.Dyn(1).B,2);
	n_x = size(in_lcsas.Dyn(1).A,1);
	n_y = size(in_lcsas.Dyn(1).C,1);
	n_w = size(in_lcsas.Dyn(1).B_w,2);
	n_v = size(in_lcsas.Dyn(1).C_v,2);

	if ~exist('debug_flag')
		debug_flag = 1; %Verbosity of Functions. Gives debugging info
	end

	cg = constr_gen();

	ops = sdpsettings('verbose',debug_flag);

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%
	if ~exist('BG')
		BG = BeliefGraph(in_lcsas,P_u,P_x0,'verbosity',debug_flag);
	end
	
	if debug_flag > 0
		figure;
		BG.plot()
	end

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
			[ Pi1{dual_var_ind} , Piu{dual_var_ind} , temp_constrs ] = cg.get_robust_reachability_constraints(	in_lcsas, possible_word, ...
																					P_x0,P_target, ...
																					Q{blang_idx},r{blang_idx}, ...
																					'P_u' , P_u);
			constraints = constraints + temp_constrs;

		end

		disp(['Created the constraints for ' num2str(blang_idx) ' belief words.' ])

	end

	%Add Block Lower Diagonal Constraints
	l_diag_constrs = cg.get_causal_constr_on_extd_gains(in_lcsas,Q);
	constraints = constraints + l_diag_constrs;

	%Insert prefix constraint
	pref_constrs = cg.create_prefix_constr_on_gains( in_lcsas , Q , r );
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
			[H0,S0,Cm0,J0,f_bar,B_w_big,C_v_big] = in_lcsas.get_mpc_matrices('word',poss_words.words{1});

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

		contr = FHAE_pb(BG,F_set,u0_set);
	end

end