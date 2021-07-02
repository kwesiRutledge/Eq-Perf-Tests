function [containment_matrix] = internal_behavior_sets2containment_mat( varargin )
	%Description:
	%
	%
	%Usage:
	%	[containment_matrix] = bg.internal_behavior_sets2containment_mat( internal_behavior_sets_in )
	%	[containment_matrix] = bg.internal_behavior_sets2containment_mat( internal_behavior_sets_in , 'verbosity' , 0)
	%	
	%

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	belief_graph_in = varargin{1};
	internal_behavior_sets_in = varargin{2};

	arg_idx = 3;
	while arg_idx <= nargin
		switch varargin{arg_idx}
			case 'verbosity'
				verbosity = varargin{arg_idx+1};
				arg_idx = arg_idx + 2;
			otherwise
				error('Unrecognized input to the function.')
		end
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	System = belief_graph_in.lcsas;

	n_sets = length(internal_behavior_sets_in);

	lcsas_in = belief_graph_in.lcsas;
	fb_method = belief_graph_in.FeedbackMethod;

	ib_dim = internal_behavior_sets_in(1).Dim;

	%Find the value of time horizon t
	lcsas_in = belief_graph_in.lcsas;
	[ n_x , n_u , n_y , n_w , n_v ] = lcsas_in.Dimensions();

	first_ibs = internal_behavior_sets_in(1);
	t = first_ibs.t;

	external_beh_dim = n_y*(t+1)+n_u*t;

	cg = constr_gen(0);

	if ~exist('verbosity')
		verbosity = 1;
	end

	ops0 = sdpsettings(	'verbose',verbosity,'cachesolvers',1,...
						'solver','gurobi', 'gurobi.BarIterLimit', 15000);

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	containment_matrix = false(length(internal_behavior_sets_in));

	for set_idx = 1:n_sets
		containment_matrix(set_idx,set_idx) = true;
	end

	for x_idx = 1:n_sets
		potential_y_idcs = [1:size(containment_matrix,2)];
		potential_y_idcs = potential_y_idcs( potential_y_idcs ~= x_idx );
		for y_idx = potential_y_idcs
			%Set up optimization
			ib_set_x = internal_behavior_sets_in(x_idx);
			ib_set_y = internal_behavior_sets_in(y_idx);

			dim_ibX = ib_set_x.Dim;
			dim_ibY = ib_set_y.Dim;

			Rt_X = [ eye(external_beh_dim), zeros(external_beh_dim, dim_ibX - external_beh_dim) ];
			Rt_Y = [ eye(external_beh_dim), zeros(external_beh_dim, dim_ibY - external_beh_dim) ];

			[dual_vars, constr] = cg.create_sadraddini_AH_inclusion_constr( ...
									zeros(external_beh_dim,1) , Rt_X , [ib_set_x.A;ib_set_x.Ae;-ib_set_x.Ae] , [ib_set_x.b;ib_set_x.be;-ib_set_x.be] , ...
									zeros(external_beh_dim,1) , Rt_Y , [ib_set_y.A;ib_set_y.Ae;-ib_set_y.Ae] , [ib_set_y.b;ib_set_y.be;-ib_set_y.be] );

			diagnostics = optimize(constr,[],ops0);

			containment_matrix(x_idx,y_idx) = (diagnostics.problem == 0);

		end
	end


end