function [ fb_law , opt_out , constraints ] = solve_reachability_problem(varargin)
%Description:
%	Synthesizes a controller that solves the reachability problem defined by the system LCSAS 
%	(including initial state set) and the target set XT.
%
%Usage:
%	[fb_law,opt_out] = lcsas.solve_reachability_problem(P_target)
%	[fb_law,opt_out,constraints] = lcsas.solve_reachability_problem(P_target)
%	[fb_law,opt_out] = lcsas.solve_reachability_problem(P_target,'fb_type','Disturbance')
%	[fb_law,opt_out] = lcsas.solve_reachability_problem(P_target,'P_u',P_u)


	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if nargin < 2
		error('Expected at least 2 inputs. Please see usage comment. (help LCSAS.rec_synthesis)')
	end

	lcsas = varargin{1};
	P_x0 = varargin{2};

	argin_idx = 3;
	while argin_idx <= nargin
		switch varargin{argin_idx}
			case 'P_u'
				P_u = varargin{argin_idx+1};
				argin_idx = argin_idx + 2;
			case 'fb_type'
				%Specifies the feedback type desired for this problem
				fb_type = varargin{argin_idx+1};
				argin_idx = argin_idx + 2;
			case 'verbosity'
				verbosity = varargin{argin_idx};
				argin_idx = argin_idx + 2;
			otherwise
				error(['Unexpected input to the function!: ' varargin{argin_idx} ])
		end
	end

	if ~exist('verbosity')
		verbosity = 0;
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	[ n , m , p , wd , vd ] = lcsas.Dimensions();

	select_m = @(t,T_r) [zeros(n,t*n), eye(n), zeros(n,(T_r-t)*n) ];

	ops = sdpsettings('verbose',verbosity);

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	%Create Belief Language through the BeliefGraph
	BG = BeliefGraph(lcsas, lcsas.L, P_u, P_x0,'verbosity',verbosity);

	cg = constr_gen();

	Q = {}; r = {}; dual_vars = {};
	constraints = [];
	for bpath_idx = 1 : length(BG.BeliefLanguage.words)
		%Get Belief Path
		bpath = BG.BeliefLanguage.words{bpath_idx};
		%Synthesize Gains and Constraints associated with this belief path
		[ Q{bpath_idx} , r{bpath_idx} , dual_vars{bpath_idx} , temp_constrs ] = BG.create_path_constraints( bpath , P_x0, P_x0 , P_u , L );
		constraints = constraints + temp_constrs;
	end

	%Add Prefix Constraints
	pref_constrs = cg.create_prefix_constr_on_gains( lcsas , BG.BeliefLanguage , F , f );
	constraints = constraints + pref_constrs;

	%% Call Optimizer %%
	optim0 = optimize(constraints,[],ops);

	%% Get Feedback Law 

	Q_set = {}; r_set = {};
	for pattern_ind = 1 : length(BG.BeliefLanguage.words)
		Q_set{pattern_ind} = value(Q{pattern_ind});
		r_set{pattern_ind} = value(r{pattern_ind});
	end

	fb_law = POB_Feedback(BG,Q_set,r_set);
	
	%% Optimization flags
	opt_out = optim0;
	opt_out.F_set = F_set;
	opt_out.f_set = f_set;


end