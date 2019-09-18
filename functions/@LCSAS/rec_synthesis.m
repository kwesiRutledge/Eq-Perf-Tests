function [ fb_law , opt_out , constraints ] = rec_synthesis(varargin)
%Description:
%	Synthesizes a controller that verifies if the system can guarantee a return to an initial set X0 after the switching behavior
%	of the input language is done.
%
%Usage:
%	[fb_law,opt_out] = lcsas.rec_synthesis(P_x0)
%	[fb_law,opt_out] = lcsas.rec_synthesis(P_x0,'P_u',P_u)
%	[fb_law,opt_out,constraints] = lcsas.rec_synthesis(P_x0,'P_u',P_u)


	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if nargin < 2
		error('Expected at least 3 inputs. Please see usage comment. (help LCSAS.rec_synthesis)')
	end

	lcsas = varargin{1};
	P_x0 = varargin{2};

	argin_idx = 3;
	while argin_idx <= nargin
		switch varargin{argin_idx}
			case 'P_u'
				P_u = varargin{argin_idx+1};
				argin_idx = argin_idx + 2;
			case 'verbosity'
				verbosity = varargin{argin_idx};
				argin_idx = argin_idx + 2;
			otherwise
				error(['Unexpected input to the function!: ' varargin{argin_idx} ])
		end
	end

	if ~exist('verbosity')
		verbosity = 1;
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	n = size(lcsas.Dyn(1).A,1);
	m = size(lcsas.Dyn(1).B,2);
	p = size(lcsas.Dyn(1).C,1);
	wd = size(lcsas.Dyn(1).B_w,2); %Assume that the size of the disturbance doesn't change for any of the included dynamics
	vd = size(lcsas.Dyn(1).C_v,2);

	select_m = @(t,T_r) [zeros(n,t*n), eye(n), zeros(n,(T_r-t)*n) ];

	ops = sdpsettings('verbose',verbosity);

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	%Create Belief Language through the BeliefGraph
	BG = BeliefGraph(lcsas, lcsas.L, P_u, P_x0,'verbosity',1);

	cg = constr_gen();

	Q = {}; r = {}; dual_vars = {};
	constraints = [];
	for bpath_idx = 1 : length(BG.BeliefLanguage.words)
		%Get Belief Path
		bpath = BG.BeliefLanguage.words{bpath_idx};
		%Synthesize Gains and Constraints associated with this belief path
		[ Q{bpath_idx} , r{bpath_idx} , dual_vars{bpath_idx} , temp_constrs ] = BG.create_path_constraints( bpath , P_x0, P_x0 , P_u);
		constraints = constraints + temp_constrs;
	end

	%Add Prefix Constraints
	pref_constrs = cg.create_prefix_constr_on_gains( lcsas , BG.BeliefLanguage , Q , r );
	constraints = constraints + pref_constrs;

	%% Call Optimizer %%
	optim0 = optimize(constraints,[],ops);

	%% Get Feedback Law 

	Q_set = {}; r_set = {};
	for pattern_ind = 1 : length(BG.BeliefLanguage.words)
		Q_set{pattern_ind} = value(Q{pattern_ind});
		r_set{pattern_ind} = value(r{pattern_ind});
	end

	fb_law = FHAE_pb(BeliefLanguage,lcsas,Q_set,r_set);
	
	%% Optimization flags
	opt_out = optim0;
	opt_out.Q_set = Q_set;
	opt_out.r_set = Q_set;


end