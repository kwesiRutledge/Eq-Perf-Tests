function [results] = observer_comparison101( varargin )
	%observer_comparison101.m
	%Description:
	%	Attempting to provide a small example showing where YALMIP's Gurobi interface appears to be wrong.

	disp(' ')
	disp('Beginning observer_comparison101.m')
	disp('Also, we are not using the linear, memory-full gains. Just offsets.')
	disp(' ')

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%
    eps0 = 10^(-5);
	TimeHorizon = 4;
	
	eta_u = 0.5;
	UPolytope = struct( ...
		'A', [1;-1] , ...
		'b', eta_u*ones(2,1) ...
		);

	eta_w = 0.25;
	WPolytope = struct( ...
		'A', [1;-1] , ...
		'b', eta_w*ones(2,1) ...
		);
	
	eta_x0 = 0.25;
	X0Polytope = struct( ...
		'A', [1;-1] , ...
		'b', eta_x0*ones(2,1) ...
		);

	x0 = -eta_x0;
	
	TargetPolytope1 = struct( ...
		'A', [1;-1] , ...
		'b', [ TimeHorizon*eta_u + eta_x0 + TimeHorizon*eta_w; -(TimeHorizon*eta_u -eta_x0 - TimeHorizon*eta_w) ] ...
		);
    TargetPolytope2 = struct( ...
		'A', [1;-1] , ...
		'b', [ (TimeHorizon-2)*eta_u + eta_x0 + TimeHorizon*eta_w; -(TimeHorizon*eta_u -eta_x0 - TimeHorizon*eta_w) ] ...
		);


	A = 1;
	B = 1;

	n_x = 1; %Dimension of x
	n_w = 1; %Dimension of w
	n_u = 1; %Dimension of u

	results.Parameters.A = A;
	results.Parameters.B = B;
	results.Parameters.n_x = n_x;
	results.Parameters.n_w = n_w;

	results.Parameters.UPolytope = UPolytope;
	results.Parameters.WPolytope = WPolytope;
	results.Parameters.X0Polytope = X0Polytope;
	results.Parameters.TargetPolytope = TargetPolytope1;

	select_state_at_T = [zeros(n_x,n_x*TimeHorizon),eye(n_x)];
	select_states_from_0_to_T = [ eye(n_x*TimeHorizon) , zeros(n_x*TimeHorizon,n_x) ];

	%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Create MPC Matrices %%
	%%%%%%%%%%%%%%%%%%%%%%%%%

	S_x0 = [];
	for t = 0:TimeHorizon
		S_x0 = [S_x0; A^t];
	end

	S_w = zeros(n_x,n_w*TimeHorizon);
	for t = 1:TimeHorizon
		S_w = [ S_w ;
				A^(t-1), S_w(end-[n_x-1:0],[1:n_w*(TimeHorizon-1)]) ];
	end

	S_u = zeros(n_x,n_u*TimeHorizon);
	for t = 1:TimeHorizon
		S_u = [ S_u ;
				A^(t-1)*B, S_u(end-[n_x-1:0],[1:n_w*(TimeHorizon-1)]) ];
	end

	results.Parameters.S_x0 = S_x0;
	results.Parameters.S_w = S_w;
	results.Parameters.S_u = S_u;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Set Up Optimization Variables %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Create Variables
	u_star1 = sdpvar(TimeHorizon*n_u,1,'full');
	u_star2 = sdpvar(TimeHorizon*n_u,1,'full');

	%%%%%%%%%%%%%%%%%%%%%%%%
	%% Set Up Constraints %%
	%%%%%%%%%%%%%%%%%%%%%%%%

	%% Constrain u_star to remain in set at all time
		
	PuT = struct( ...
		'A', kron(eye(TimeHorizon),UPolytope.A), ...
		'b', repmat(UPolytope.b,[TimeHorizon,1]) ...  
		);
	
	input_bounds_constraint = [ PuT.A*u_star1 <= PuT.b ] + [ PuT.A*u_star2 <= PuT.b ];

	%% Overlapping W's Exist

	MatchingSet_Ws = struct( ...
		'A', [	-eye(n_x*(TimeHorizon+1)), ( S_w ) * [ eye(n_w*TimeHorizon) , zeros(n_w*TimeHorizon) ] ;
				-eye(n_x*(TimeHorizon+1)), ( S_w ) * [ zeros(n_w*TimeHorizon) , eye(n_w*TimeHorizon)] ] , ...
		'b', [ 	S_x0 * x0 + S_u * u_star1 ;
				S_x0 * x0 + S_u * u_star2 ] ...
		)

	% show that there exists a disturbance which makes this set nonempty
% 	eta = sdpvar(size(MatchingSet_Ws.A,2),1,'full');
% 	nonempty_constraint = [ MatchingSet_Ws.A * eta <= MatchingSet_Ws.b ];

    y = sdpvar(size(MatchingSet_Ws.A,1),1,'full');
    empty_constraint = [MatchingSet_Ws.A'*y == 0] + [ MatchingSet_Ws.A'*y <= -eps0 ] + [y >= 0]
    
   	matching_input_constraint = [ u_star1([1:n_x*(TimeHorizon-3)]) == u_star2([1:n_x*(TimeHorizon-3)]) ];

	%% Create Robust Reachability Constraints

	PwT = struct( ...
		'A' , kron(eye(TimeHorizon),WPolytope.A) , ...
		'b' , repmat(WPolytope.b,[TimeHorizon,1]) );

	ReachableSetWs_1 = struct( ...
		'A', TargetPolytope1.A * select_state_at_T* S_w, ...
		'b', TargetPolytope1.b - TargetPolytope1.A * select_state_at_T * ( S_x0 * x0 + S_u * u_star1 ) ...
		);

	ReachableSetWs_2 = struct( ...
		'A', TargetPolytope2.A * select_state_at_T* S_w, ...
		'b', TargetPolytope2.b - TargetPolytope2.A * select_state_at_T * ( S_x0 * x0 + S_u * u_star2 ) ...
		);

	Lambda_Reach1 = sdpvar(size(ReachableSetWs_1.A,1),size(PwT.A,1),'full');
	Lambda_Reach2 = sdpvar(size(ReachableSetWs_2.A,1),size(PwT.A,1),'full');
	reachability_constraint = ...
		[Lambda_Reach1 >= 0] + [ Lambda_Reach1 * PwT.A == ReachableSetWs_1.A ] + [ Lambda_Reach1 * PwT.b <= ReachableSetWs_1.b ] + ...
		[Lambda_Reach2 >= 0] + [ Lambda_Reach2 * PwT.A == ReachableSetWs_2.A ] + [ Lambda_Reach2 * PwT.b <= ReachableSetWs_2.b ];

	%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Optimize With Gurobi %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%

	optimization_constraints = 	reachability_constraint + input_bounds_constraint +  empty_constraint;

	ops = sdpsettings('verbose',1,'debug',1);
	ops = sdpsettings(ops,'solver','gurobi');

	optim0 = optimize(optimization_constraints,[],ops);

	results.GurobiSolution.Output = optim0;
	results.GurobiSolution.u_star1 = value(u_star1);
	results.GurobiSolution.u_star2 = value(u_star2);

	% Check Solution
	disp('Checking Gurobi Solution')
	disp('Is final state in the target region for an arbitrary bad disturbance trajectory w?')
	w_bad1 = -eta_w*ones(TimeHorizon*n_w,1);
	w_bad2 = eta_w*ones(TimeHorizon*n_w,1);

	disp(all(TargetPolytope1.A* select_state_at_T * ( S_x0 * x0 + S_u * value(u_star1) + S_w * w_bad1 ) <= TargetPolytope1.b))
	disp(all(TargetPolytope1.A* select_state_at_T * ( S_x0 * x0 + S_u * value(u_star1) + S_w * w_bad2 ) <= TargetPolytope1.b))

	disp('Plotting...')
	PwT_mpt3 = Polyhedron('A',PwT.A,'b',PwT.b);
	same_axis = [-1,4,-1,1];

	figure(1)
	subplot(2,1,1)
	plot( select_state_at_T*(S_x0*x0+S_u*value(u_star1)+S_w*PwT_mpt3) )
	axis(same_axis)
	title('Reachable Set at T (Gurobi)')

	subplot(2,1,2)
	plot( Polyhedron('A',TargetPolytope1.A,'b',TargetPolytope1.b) )
	axis(same_axis)
	title('Target Set')


	disp(' ')


	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Optimize With fmincon %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%

	ops = sdpsettings('verbose',1,'debug',1);
	ops = sdpsettings(ops,'solver','fmincon');

	optim1 = optimize(optimization_constraints,[],ops);

	results.FminconSolution.Output = optim1;
	results.FminconSolution.u_star1 = value(u_star1);
	results.FminconSolution.u_star2 = value(u_star2);

	disp('Checking fmincon Solution')
	disp('Is final state in the target region for an arbitrary bad disturbance trajectory w?')
	w_bad = -eta_w*ones(TimeHorizon*n_w,1);

	disp(all(TargetPolytope1.A* select_state_at_T * ( S_x0 * x0 + S_u * value(u_star1) + S_w * w_bad ) <= TargetPolytope1.b))
	disp(all(TargetPolytope1.A* select_state_at_T * ( S_x0 * x0 + S_u * value(u_star1) + S_w * w_bad2 ) <= TargetPolytope1.b))

	disp('Numerical issues cause the result to be a bit fishy.')

	disp('Plotting...')
	PwT_mpt3 = Polyhedron('A',PwT.A,'b',PwT.b);
	same_axis = [-1,4,-1,1];

	figure(2)
	subplot(2,1,1)
	plot( select_state_at_T*(S_x0*x0+S_u*value(u_star1)+S_w*PwT_mpt3) )
	axis(same_axis)
	title('Reachable Set at T (fmincon)')

	subplot(2,1,2)
	plot( Polyhedron('A',TargetPolytope1.A,'b',TargetPolytope1.b) )
	axis(same_axis)
	title('Target Set')

end