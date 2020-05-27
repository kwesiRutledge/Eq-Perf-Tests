function [containment_matrix] = internal_behavior_sets2containment_mat( belief_graph_in , internal_behavior_sets_in )
	%Description:
	%
	%
	%Usage:
	%	[containment_matrix] = bg.internal_behavior_sets2containment_mat( t_in , internal_behavior_sets_in )
	%


	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	n_sets = length(internal_behavior_sets_in);

	lcsas_in = belief_graph_in.lcsas;
	fb_method = belief_graph_in.FeedbackMethod;

	ib_dim = internal_behavior_sets_in(1).Dim;

	%Find the value of time horizon t
	lcsas_in = belief_graph_in.lcsas;
	n_x = size(lcsas_in.Dyn(1).A,1);
	n_u = size(lcsas_in.Dyn(1).B,2);
	n_w = size(lcsas_in.Dyn(1).B_w,2);
	if strcmp(fb_method,'output')
		n_y = size(lcsas_in.Dyn(1).C,1);
		n_v = size(lcsas_in.Dyn(1).C_v,2);

		t = (ib_dim-n_y-n_v-2*n_x)/(n_y+n_u+n_w+n_v+n_x);

	elseif strcmp(fb_method,'state')
		t = (ib_dim - 2*n_x)/(n_x+n_u+n_w);
	else
		error('Unrecognized feedback method.')
	end

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