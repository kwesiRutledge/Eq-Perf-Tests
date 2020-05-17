function [results] = observer_comparison76( varargin )
	%observer_comparison76.m
	%Description:
	%	Comparing the method for detecting projection inclusion that I have against sadra's condition.
	%

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	% if nargin >= 1
	% 	load_data_flag = varargin{1};
	% end

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

	L1 = Language([1,2,1,2],[3,4,3,4],[5,1,1,1]);
	T = 4;

	%Create a simple Language Constrainted Switching System
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

	f1 = eta_w*[0;1];
	f2 = eta_w*[1;0];
	f3 = -f1;
	f4 = -f2;

	aff_dyn_list = [	Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f3,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f4,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f1,C1,Pw1,Pv2) ];

	lcsas0 = LCSAS( aff_dyn_list , L1 , 'X0' , P_x0 );

	results.lcsas = lcsas0;

	%% Debugging Variables

	verbosity = 1;
	update_freq = 50;

	results.parameters.dim = size(A1,1);

	%%%%%%%%%%%
	%% Tests %%
	%%%%%%%%%%%

	disp(' ')
	disp(['Beginning ' test_name '.' ])
	disp('The objective of this test is to create a test-bed for creating BeliefGraphs/')
	disp(' ')

	%% Get System %%
	lcsas0.X0 = P_x0;
	bg0 = BeliefGraph( lcsas0 , P_u , P_x0 , 'return_empty' , true );

	%% Create Initial Nodes %%
	N0 = bg0.get_initial_beliefnodes()
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
			temp_post = bg0.post(c_node,P_u,P_x0,'debug',verbosity);
								
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

	%Create Belief Language
	%bg0.BeliefLanguage = bg0.get_belief_language();


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Save The Larger Variables to a Data File %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%save([ save_file_name '.mat'])

end