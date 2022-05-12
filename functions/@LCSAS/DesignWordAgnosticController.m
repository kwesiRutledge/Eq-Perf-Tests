function [ controller , info ] = DesignWordAgnosticController( varargin )
	%Description
	%	Designs a controller for the system lcsas0 to reach the target P_target with:
	%	- input constraints
	%	- reaches the target
	%
	%Usage:
	%	[ contr , info ] = lcsas.DesignWordAgnosticController(P_target)

	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	[ System , P_target , settings ] = ip_DesignWordAgnosticController(varargin{:});

	%% Constants %%
	%%%%%%%%%%%%%%%

	[ n_x , n_u , n_y , n_w , n_v ] = System.Dimensions();
	L = System.L;

	X0 = System.X0; U = System.U;

	cg = constr_gen(0);

	T = length(L.words{1});

	%% Algorithm %%
	%%%%%%%%%%%%%%%

	% Define Gains
	Q = sdpvar(n_u*T,n_x*(T+1),'full');
	r = sdpvar(n_u*T,1);

	% Create Input Sequence Polytope
	U_T = 1;
	for t = 1:T 
		U_T = U_T * U;
	end

	% Create Constraints
	constraints = []
	for word_index = [1:L.cardinality()]

		% Get Word
		temp_word = L.words{word_index};

		% Create Reachability Constraint
		[ S_w , S_u , S_C , S_x0 , S_k , S_Bw , S_Cv ] = System.get_mpc_matrices('word',temp_word); %Get MPC Matrices

		I_r = eye(size(S_u,1),size(Q,2));

		G = [ (I_r+S_u*Q) * S_w , (I_r+S_u*Q) * S_x0 ];
		selectFinalState = [ zeros(n_x,n_x*T) , eye(n_x) ];

		H_T = P_target.A * selectFinalState * G;
		h_T = P_target.b - P_target.A * selectFinalState * ( S_u * r);

		P_eta = System.ToWSequence(temp_word) * X0;

		[ dual_vars , reachability_constraints ] = cg.get_H_polyt_inclusion_constr( P_eta.A , P_eta.b , H_T , h_T );
		constraints = constraints + reachability_constraints;

		% Create Input Bound Constraint

		G_u = [ Q*S_w , Q*S_x0 ];
		H_T = U_T.A * G_u;
		h_T = U_T.b - U_T.A * (r);

		[ dual_vars , input_bound_constraints ] = cg.get_H_polyt_inclusion_constr( P_eta.A , P_eta.b , H_T , h_T );
		constraints = constraints + input_bound_constraints;

	end

	causal_gain_constraints = System.create_lower_diagonal_constraint_on_gains( { Q } , 'Disturbance' );
	constraints = constraints + causal_gain_constraints;

	%% Optimize
	ops = sdpsettings('verbose',settings.verbosity,'debug',1);
	ops = sdpsettings(ops,'solver','gurobi');

	optim0 = optimize(constraints,[],ops);

	info.yalmiptime = optim0.yalmiptime;
	info.solvertime = optim0.solvertime; 

	%% Create Outputs
	Q_star = value(Q); info.Q_star = Q_star;
	r_star = value(r); info.r_star = r_star;

	controller.K = eye(size(Q_star,1),size(S_u,2)) + Q_star * S_u;
	controller.u0 = ( eye(size(Q_star,1),size(S_u,2)) + Q_star * S_u )^(-1) * r_star;

end

function [ System , P_target , settings ] = ip_DesignWordAgnosticController(varargin)
	%Description:
	%	Processes the inputs given to DesignWordAgnosticController. May produce an error if some things are not
	%	defined or provided.
	%Usage:
	%	[ System , P_target , settings ] = ip_DesignWordAgnosticController(varargin{:})

	%% Get First two items

	System = varargin{1};
	P_target = varargin{2};

	%% Check System to make sure that X0 and U sets are defined.
	System.check('X0','U')

	%% Create Default Settings
	settings = struct( ...
		'FeedbackType','Disturbance', ...
		'verbosity', 1 ...
		);

end