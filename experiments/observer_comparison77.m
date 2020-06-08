function [results] = observer_comparison77( varargin )
	%observer_comparison76.m
	%Description:
	%	Comparing the method for detecting projection inclusion that I have against sadra's condition.
	%

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if nargin >= 1
		experimental_post_flag = varargin{1};
	end

	if nargin >= 2
		system_name = varargin{2};
	end

	if nargin >= 3
		verbosity = varargin{3};
	end

	% if nargin >= 3
	% 	c_sq.dim_x = varargin{2};
	% 	c_sq.dim_y = varargin{3};
	% end

	% if nargin >= 4
	% 	verbosity = varargin{4};
	% end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	test_name = 'observer_comparison77';
	save_file_name = [ 'results/' test_name '_results.mat'];

	%% Create System

	if ~exist('system_name')
		system_name = 'simple_2d_integrator';
	end

	%Create a simple Language Constrainted Switching System
	switch system_name
		case 'simple_2d_integrator'
			L1 = Language([1,2,1,2],[3,4,3,4],[5,1,1,1]);
			T = length(L1.words{1});

			A1 = [0,1;0.1,-0.05];
			B1 = [0;1];
			C1 = [1,0];
			
			n_x = size(A1,1);
			n_u = size(B1,2);
			n_y = size(C1,1);

			eta_v = 0.1; eta_w = 0.2;
			Pv1 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));
			Pw1 = Polyhedron('lb',-eta_w*ones(1,n_x) ,'ub',eta_w*ones(1,n_x));
			eta_v = 0.3;
			Pv2 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));

			eta_u = 0; eta_x0 = 0.3;
			P_u = Polyhedron('lb',-eta_u*ones(1,n_u) ,'ub',eta_u*ones(1,n_u));
			P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));

			f1 = eta_w*[zeros(n_x-1,1);1];
			f2 = eta_w*[1;zeros(n_x-1,1)];
			f3 = -f1;
			f4 = -f2;

			aff_dyn_list = [	Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1), ...
								Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1), ...
								Aff_Dyn(A1,B1,f3,C1,Pw1,Pv1), ...
								Aff_Dyn(A1,B1,f4,C1,Pw1,Pv1), ...
								Aff_Dyn(A1,B1,f1,C1,Pw1,Pv2) ];

			lcsas0 = LCSAS( aff_dyn_list , L1 , 'X0' , P_x0 );

		case 'rand_4d_sys'
			L1 = Language([1,2,1,2],[3,4,3,4],[5,1,1,1]);
			T = length(L1.words{1});

			A1 = randn(4);
			B1 = [zeros(2);eye(2)];
			C1 = [1,0,0,0;0,0,1,0];
			
			n_x = size(A1,1);
			n_u = size(B1,2);
			n_y = size(C1,1);

			eta_v = 0.1; eta_w = 0.2;
			Pv1 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));
			Pw1 = Polyhedron('lb',-eta_w*ones(1,n_x) ,'ub',eta_w*ones(1,n_x));
			eta_v = 0.3;
			Pv2 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));

			eta_u = 0; eta_x0 = 0.3;
			P_u = Polyhedron('lb',-eta_u*ones(1,n_u) ,'ub',eta_u*ones(1,n_u));
			P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));

			f1 = eta_w*[zeros(n_x-1,1);1];
			f2 = eta_w*[1;zeros(n_x-1,1)];
			f3 = -f1;
			f4 = -f2;

			aff_dyn_list = [	Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1), ...
								Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1), ...
								Aff_Dyn(A1,B1,f3,C1,Pw1,Pv1), ...
								Aff_Dyn(A1,B1,f4,C1,Pw1,Pv1), ...
								Aff_Dyn(A1,B1,f1,C1,Pw1,Pv2) ];

			lcsas0 = LCSAS( aff_dyn_list , L1 , 'X0' , P_x0 );

		case 'oc60_L6'
			L = Language([1,2,2,1,1,1],[1,2,1,1,1,1],[1,1,2,1,1,1]);
			T = length(L.words{1});

			c_sq.dim_x = 2;
			c_sq.dim_y = 1;
			dt = 0.1;

			eta_w = 0.35; eta_v = 0.2;

			in_sys = get_consensus_dyn(c_sq.dim_x,c_sq.dim_y,dt,'L',L,'disturb_params', eta_w , eta_v );

			n_x = size(in_sys.Dyn(1).A,1); n_u = size(in_sys.Dyn(1).B,2);
			eta_x0 = 0.3; eta_u = 50;
			P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));
			P_u = Polyhedron('lb',-eta_u*ones(1,n_u),'ub',eta_u*ones(1,n_u));

			in_sys.X0 = P_x0;

			lcsas0 = in_sys;

		case 'oc60_L4'
			L = Language([1,2,2,1],[1,2,1,1],[1,1,2,1]);
			T = length(L.words{1});

			c_sq.dim_x = 2;
			c_sq.dim_y = 1;
			dt = 0.1;

			eta_w = 0.35; eta_v = 0.2;

			in_sys = get_consensus_dyn(c_sq.dim_x,c_sq.dim_y,dt,'L',L,'disturb_params', eta_w , eta_v );

			n_x = size(in_sys.Dyn(1).A,1); n_u = size(in_sys.Dyn(1).B,2);
			eta_x0 = 0.3; eta_u = 50;
			P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));
			P_u = Polyhedron('lb',-eta_u*ones(1,n_u),'ub',eta_u*ones(1,n_u));

			in_sys.X0 = P_x0;

			lcsas0 = in_sys;

		case 'compass_L4'
			%Define constants from Sadra's file
			m=2; m_h=2; a=0.5; b=0.5;
			n = 4;
			dt=0.04; 
			T = 10; %T=45;
			g = 10;
			l = a+b;

			M = [ 	[(m_h+m)*(l^2) + m*a^2, m*l*b];
					[ m*l*b,m*b^2] ];
			Minv = pinv(M);
			tau_g = g*[ [m_h*l+m*a+m*l,0] ;
						[0,-m*b] ];
			aff_dyn_list = [];
			for t = 0:T-1
				A_t = eye(n);
				A_t([1:2],[3:4]) = eye(2)*dt;
				A_t([3:end],[1:2]) = Minv*tau_g*dt;
				% A[t][0,2],S.A[t][1,3]=dt,dt
				%S.A[t][2:,0:2]=np.dot(Minv,tau_g)*dt
				Bq = (Minv * ones(2,1))*dt;
				B_t = [0;0;Bq];
				% Bq = np.dot(Minv,np.array([1,1]))*dt
				% S.B[t]=np.array([0,0,Bq[0],Bq[1]]).reshape(4,1)
				
				% S.B[t]=np.array([0,0,1,1]).reshape(4,1)
				% S.B[t]=np.array([[0,0,1,0],[0,0,0,1]]).T
				% S.B[t]=np.array([0,0,1,1]).reshape(4,1)
				% S.C[t]=np.array([[1,0,0,0],[0,1,0,0]]).reshape(2,4)
			    C_t = [1,1,0,0];
			    n_y = size(C_t,1);
			    %S.C[t]=np.array([[1,1,0,0]]).reshape(1,4)
			    D = eye(n);
			    %S.D[t]=np.eye(n)
			    d = zeros(n,1);
			  	%S.d[t]=np.zeros((n,1))
			  	eta_w = 10^(-3);
			  	P_w = Polyhedron('lb',-eta_w*ones(1,n),'ub',eta_w*ones(1,n));
			    %S.W[t]=zonotope(np.ones((n,1))*0,np.ones((n,1))*0)
				%S.V[t]=zonotope(np.zeros((2,1)),np.ones((2,1))*0.0000)
			    eta_v = 10^(-4);
			    P_v = Polyhedron('lb',-eta_v*ones(1,n_y),'ub',eta_v*ones(1,n_y));
			    %S.V[t]=zonotope(np.zeros((1,1)),np.eye(1)*0.00)
			    aff_dyn_list = [aff_dyn_list,Aff_Dyn( A_t,B_t,zeros(n,1),C_t,P_w,P_v )];
			end

			L = Language([1,2,3,4],[1,2,2,3],[1,2,3,3]);
			T = length(L.words{1});

			%Define P_x0 and P_u
			n_x = size(A_t,1);
			n_u = size(B_t,2);
			eta_x0 = 10^(-4); eta_u = 50;
			P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));
			P_u = Polyhedron('lb',-eta_u*ones(1,n_u),'ub',eta_u*ones(1,n_u));

			lcsas0 = LCSAS( aff_dyn_list , L , 'X0' , P_x0 );

		otherwise
			error('Unexpected system name.')
	end
	

	results.lcsas = lcsas0;

	if ~exist('experimental_post_flag')
		experimental_post_flag = 0;
	end

	%% Debugging Variables
	if ~exist('verbosity')
		verbosity = 1;
	end

	update_freq = 50;

	results.parameters.dim = n_x;

	%%%%%%%%%%%
	%% Tests %%
	%%%%%%%%%%%

	disp(' ')
	disp(['Beginning ' test_name '.' ])
	disp('The objective of this test is to create a test-bed for creating BeliefGraphs/')
	disp(' ')
	disp('Parameters:')
	disp(['experimental_post_flag = ' num2str(experimental_post_flag)])
	disp(['system_name = ' system_name ])
	disp(['verbosity = ' num2str(verbosity)])
	disp(' ')

	mptopt('lpsolver','mosek');

	%% Get System %%
	bg_construction_start = tic;

	lcsas0.X0 = P_x0;
	if (experimental_post_flag == 1) | (experimental_post_flag == 4)
		bg0 = BeliefGraph( lcsas0 , P_u , P_x0 , 'return_empty' , true , ...
											 	'use_proj_flag', true , ...
											 	'ConsistencySetVersion' , 2 );
	else	
		bg0 = BeliefGraph( lcsas0 , P_u , P_x0 , 'return_empty' , true , ...
												 'use_proj_flag', true );
	end

	%% Create Initial Nodes %%
	N0 = bg0.get_initial_beliefnodes('OverrideProjFlag',false);
	results.N0 = N0;

	nodes_at_time_tm1 = N0;
	bg0.N = N0;

	for tau = 1:T-1
		%Each belief will be indexed by a time. (i.e. I hold X belieft at time t)

		if verbosity > 0
			disp(['- tau = ' num2str(tau) ])
		end

		for node_idx = 1:length(nodes_at_time_tm1) %Iterate through all nodes that are stored in the nodes_at_time_tm1 array
			%Current node
			c_node = nodes_at_time_tm1(node_idx);

			%Calculate the ancestors of this Belief Node
			switch experimental_post_flag
				case 0
					temp_post = bg0.post(c_node,P_u,P_x0,'debug',verbosity);
				case 1
					temp_post = bg0.post_experimental(c_node,P_u,P_x0,'debug',verbosity);
				case 2
					temp_post = bg0.post_experimental2(c_node,P_u,P_x0,'debug',verbosity);
				case 3
					temp_post = bg0.post_accel(c_node,P_u,P_x0,'debug',verbosity);
				case 4
					temp_post = bg0.post_experimental2_accel(c_node,P_u,P_x0,'debug',verbosity);
				case 5
					temp_post = bg0.post_experimental_accel(c_node,P_u,P_x0,'debug',verbosity);
				case 6
					temp_post = bg0.post_experimental_accel(c_node,P_u,P_x0,'debug',verbosity,'PerformBoundingBoxCheck',true);
				otherwise
					error(['Unexpected value of experimental_post_flag: ' num2str(experimental_post_flag) ])
			end
								
			%Add the ancestors to the BeliefGraph's set of nodes if they don't already exist in the set.
			for node_idx = 1:length(temp_post)
				if bg0.find_node_idx(temp_post(node_idx)) == -1
					bg0.N = [bg0.N,temp_post(node_idx)];
				end
			end

			% Create edges corresponding to the ancestor-descendant relationship
			for edge_idx = 1:length(temp_post)
				%Create edge using this new node and add it to the edges list
				temp_edge = [bg0.find_node_idx(c_node), ...
							 bg0.find_node_idx(temp_post(edge_idx))];
				bg0.E = [bg0.E;temp_edge];
			end

		end

		%Create next level of the tree
		nodes_at_time_tm1 = bg0.get_all_nodes_at_time(tau);
		if verbosity >= 1
			disp(['There are ' num2str(length(nodes_at_time_tm1)) ' nodes at time tau = ' num2str(tau) '.' ])
		end
	end

	results.BeliefGraph = bg0;
	results.BGConstructionTime = toc(bg_construction_start);

	%Create Belief Language
	%bg0.BeliefLanguage = bg0.get_belief_language();


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Save The Larger Variables to a Data File %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%save([ save_file_name '.mat'])

end