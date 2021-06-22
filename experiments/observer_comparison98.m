function [results] = observer_comparison98( varargin )
	%observer_comparison98.m
	%Description:
	%	Attempting to implement a solution to subproblem 3 of the LCSLS wiki pages
	%	using an enumeration based version of the approach given in observer_comparison86.m.
	%	This method also restricts itself by not using linear gains.

	disp(' ')
	disp('Beginning observer_comparison98.m')
	disp('Directly running the enumeration approach for similar rotation system')
	disp('with a modified yalmip2gurobi.m file.')
	disp('Also, we are not using the linear, memory-full gains.')
	disp(' ')

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if nargin >= 1
		TimeHorizon = varargin{1};
	end

	if nargin >= 2
		solution_approach = varargin{2};
	end

	if nargin >= 3
		system_option_val = varargin{3};
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	% Save To Parameters File

	eps0 = 10^(-5);
	bigM = 10^(5);

	temp_K_bound = 10^6;

	%P_target = Polyhedron('lb',-0.4*ones(1,2),'ub',0.4*ones(1,2));

	if ~exist('TimeHorizon')
		TimeHorizon = 4;
	end

	if ~exist('solution_approach')
		solution_approach = 'Bilinear (Gurobi)';
	end

	if ~exist('system_option_val')
		system_option_val = 1;
	end

	RelaxationOrder = 1;

	%pool = parpool(4);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Convert LIPM Models to Aff_Dyn Objects %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	switch system_option_val
		case 1
			[ lcsas0 , TimeHorizon , Pu , Pw , x0 , Px0 , P_target ] = get_similar_rotation_lcsas('TimeHorizon',TimeHorizon);
		case 2
			[ lcsas0 , TimeHorizon , Pu , Pw , x0 , Px0 , P_target ] = get_opposing_rotation_lcsas('TimeHorizon',TimeHorizon);
		case 3
			[ lcsas0 , x0 , TimeHorizon , P_target ] = get_differently_loaded_drone_lcsas('TimeHorizon',TimeHorizon);
			Pu = lcsas0.U;
			Px0 = lcsas0.X0;
		otherwise
			error(['Expected a value between 1 and 3. Received ' num2str(system_option_val) '.' ])
	end

	results.Parameters.LCSAS = lcsas0;
	results.Parameters.P_target = P_target;

	disp('Created System Object.')

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Create All Possible Belief Sequences %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	disp('Creating All Possible Belief Sequences.')

	pathConstructionStartTime = tic();

	% This is really just the BeliefGraph algorithm but 
	% without computing the safe set.

	LK_t0 = [ lcsas0.L ];

	[ LK_sequences , LK_history ] = LK_t0.create_belief_sequences_of_length(TimeHorizon);
	num_knowl_sequences = size(LK_sequences,2);

	disp(['- Created all sequences of length ' num2str(TimeHorizon) '.' ])

	results.LK = LK_sequences;

	% Extract Disturbances Associated with Each Sequence
	LK_w = {}; 
	% for LK_sequence_index = 1:num_knowl_sequences
	% 	temp_LK_sequence = LK_sequences(:,LK_sequence_index);

	% 	[ LK_w{LK_sequence_index} ] = lcsas0.find_hypothesis_generating_disturbances( temp_LK_sequence([2:end]) );
	% 	LK_w{LK_sequence_index}.outerApprox;

	% end
	% disp(['- Collected all disturbances for each of the ' num2str(num_knowl_sequences) ' sequences.' ])

	results.LK_Sequences = LK_sequences;
	results.LK_w = LK_w;

	[ possibleSubsetsOfPaths , possible_choices , choices_as_binary_flags ] = lcsas0.get_feasible_combinations_of_beliefs( LK_sequences([2:end],:) , 'verbosity' , 1 );

	results.MatchingBehaviorInstances = choices_as_binary_flags;
	results.possibleSubsetsOfPaths = possibleSubsetsOfPaths;

	results.PathConstructionTime = toc(pathConstructionStartTime);

	disp('- Found all feasible combinations of beliefs.')

	%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Set Up Optimization %%
	%%%%%%%%%%%%%%%%%%%%%%%%%

	[ n_x , n_u , n_y , n_w , n_v ] = lcsas0.Dimensions();

	% active_gains = binvar(TimeHorizon,size(LK{TimeHorizon},2),'full');

	% matching_behavior = {};
	% for t = 1:TimeHorizon
	% 	matching_behavior{t} = binvar(size(LK{t},2),1,'full');
	% end

	OverallRuntimeStart = tic();

	for choice_index = 1:length(choices_as_binary_flags)

		%Start timer
		choice_index_start_time = tic;

		% Variables

		% Create Optimization Variables
		K = {}; k = {};
		for knowl_seq_index = 1:num_knowl_sequences
			K{knowl_seq_index} = zeros(TimeHorizon*n_u,TimeHorizon*n_w);%sdpvar(TimeHorizon*n_u,(TimeHorizon)*n_w,'full');
			k{knowl_seq_index} = sdpvar(TimeHorizon*n_u,1,'full');
		end

		matching_behavior = choices_as_binary_flags{choice_index};

		% contains_u = binvar(TimeHorizon,size(LK{TimeHorizon},2),'full');

		% Create Sets

		PuT = 1;
		for t = 1:TimeHorizon
			PuT = PuT * Pu;
		end

		PwT = {};
		for word_index = 1:lcsas0.L.cardinality()
			temp_word = lcsas0.L.words{word_index};
			PwT{word_index} = 1;
			for t = 1:TimeHorizon
				PwT{word_index} = PwT{word_index} * lcsas0.Dyn( temp_word(t) ).P_w;
			end
		end

		[ S_w , S_u , ~ , J , f_bar ] = lcsas0.get_mpc_matrices('All Words');

		% Constraints

		cg = constr_gen(0);

		%% Input Bound Constraints

		input_bounds_constraint = cg.get_input_bound_constraint_on( lcsas0 , K , k , 'fb_type' , 'state-disturbance' ); 

		%% Matching Behavior Constraints

		matching_behavior_constraint = [];
		dummy_var_bound_constraints = []; temp_dummy_w = {}; temp_dummy_y = {};
		for knowl_seq_index = 1:num_knowl_sequences
			
			knowl_seq = LK_sequences(:,knowl_seq_index);
			temp_last_lang = knowl_seq(end);
			tll_card = temp_last_lang.cardinality();

			% Constrain all controllers K{} with this prefix
			%relevant_gain_indices = find_all_knowledge_sequences_with_prefix(LK{TimeHorizon},knowl_seq);
			[ H_PhiI , h_PhiI ] = consistent_set_matrices( lcsas0 , TimeHorizon-1 , knowl_seq(end) , Pu , Px0 );

			%Get Closed Loop Consistent Internal Behavior Set Matrices (Function of K)
			[ H_cl , h_cl ] = lcsas0.get_closed_loop_consistent_internal_behavior_set_matrices( ...
				H_PhiI , h_PhiI , ...
				x0 , ...
				K{knowl_seq_index} , k{knowl_seq_index} , knowl_seq );

			 % Create Constraints
			 if matching_behavior(knowl_seq_index)
			 % 	Pw_prime = 1;
			 % 	for Pw_index = 1:size(LK_w{knowl_seq_index})
			 % 		Pw_prime = Pw_prime * LK_w{knowl_seq_index}(Pw_index);
			 % 	end

				% %Apply the existence of w condition.
				% temp_dummy_w{end+1} = sdpvar(size(H_cl,2),1,'full');
				% matching_behavior_constraint = matching_behavior_constraint + [ H_cl * temp_dummy_w{end} <= h_cl ] + [ Pw_prime.A * temp_dummy_w{end} <= Pw_prime.b ];

				%Guarantee that no other controllers paths are touched by this controller
				[ constrain_other_beliefs_tobe_empty , y_arr ] = cg.no_other_belief_sequences_reached( ...
					LK_sequences , matching_behavior , ...
					lcsas0 , ...
					K{knowl_seq_index} , k{knowl_seq_index} , ...
					x0 );

			end

		end

		%% Create Causal Detection Constraints %%
		causal_detection_constraints = cg.get_belief_prefix_gain_constraints( lcsas0 , K , k , LK_sequences , 'Ignore K');

		%% Create Robust Reachability Constraints

		guaranteed_reachability_constraint = [];
		incl_constraints = {};
		reachability_dual_vars = {};
		reachability_dual_var_constraints = [];
        
		for knowl_seq_index = 1:num_knowl_sequences
            
            H = {}; h = {};
            nonK_prefactor = {}; K_prefactor = {};
            h_independent_factor = {}; h_dependent_factor = {};
			
			knowl_seq = LK_sequences(:,knowl_seq_index);
			temp_last_lang = knowl_seq(end);

			
			%HACK
			%	Extend knowledge sequence.
			extended_knowl_seqs = [];
			for word_index = 1:temp_last_lang.cardinality()
				extended_knowl_seq = [ knowl_seq ; Language( temp_last_lang.words{word_index} ) ];
				extended_last_lang = extended_knowl_seq(end);

				%Get Closed Loop Consistent Internal Behavior Set Matrices (Function of K)
				[ H_PhiI , h_PhiI ] = consistent_set_matrices( lcsas0 , TimeHorizon , extended_knowl_seq(end) , Pu , Px0 );

				[ H_cl , h_cl ] = lcsas0.get_closed_loop_consistent_internal_behavior_set_matrices( ...
					H_PhiI , h_PhiI , ...
					x0 , ...
					K{knowl_seq_index} , k{knowl_seq_index} , extended_knowl_seq );

				switch extended_last_lang.cardinality()
				 	case 1

				 		[~,word_id] = lcsas0.L.contains(extended_last_lang.words{1});

						nonK_prefactor{1} = [ 	S_w{word_id} ;
							zeros(n_u*TimeHorizon,n_w*TimeHorizon) ;
							eye(n_w*TimeHorizon) ;
							zeros(n_x,n_w*TimeHorizon) ];

						K_prefactor{1} = [ S_u{word_id} ;
										eye(n_u*TimeHorizon);
										zeros(n_w*TimeHorizon,n_u*TimeHorizon);
										zeros(n_x,n_u*TimeHorizon) ];

						H{1} = H_PhiI * ( nonK_prefactor{1} );

						h_independent_factor{1} = h_PhiI - ...
							H_PhiI * [ J{word_id}*x0 + S_w{word_id}*f_bar{word_id} ;
									zeros(n_u*TimeHorizon,1) ;
									zeros(n_w*TimeHorizon,1) ; 
									x0 ];

						h_dependent_factor{1} = - ...
							H_PhiI * [ S_u{word_id} ;
									eye(n_u*TimeHorizon) ;
									zeros(n_w*TimeHorizon,n_u*TimeHorizon) ; 
									zeros(n_x,n_u*TimeHorizon) ];

						h{1} = h_independent_factor{1} + h_dependent_factor{1}*k{knowl_seq_index};

						[ reachability_dual_vars{end+1} , incl_constraints{end+1} ] = cg.get_H_polyt_inclusion_constr( ...
							H{1} , h{1} , ...
							P_target.A*[ zeros(n_x,n_x*TimeHorizon), eye(n_x) ]*( S_w{word_id} ) , ...
							P_target.b - P_target.A*[ zeros(n_x,n_x*TimeHorizon), eye(n_x) ]*( J{word_id}*x0 + S_u{word_id}*k{knowl_seq_index} + S_w{word_id}*f_bar{word_id} ) );

						% guaranteed_reachability_constraint = guaranteed_reachability_constraint + ...
						% 	iff( matching_behavior{TimeHorizon-1}(knowl_seq_index) == 1 , incl_constraints ) ;
						
				 	case 2
				 		for tll_index = 1:temp_last_lang.cardinality()
					 		[~,word_id] = lcsas0.L.contains(temp_last_lang.words{tll_index});

					 		nonK_prefactor{tll_index} = ...
					 			[ zeros(n_x*(TimeHorizon+1),(tll_index-1)*n_w*TimeHorizon), S_w{word_id}, zeros(n_x*(TimeHorizon+1),(tll_card-tll_index)*n_w*TimeHorizon) ;
									zeros(n_u*TimeHorizon,tll_card*n_w*TimeHorizon) ;
									eye(tll_card*n_w*TimeHorizon) ;
									zeros(tll_card*n_x,tll_card*n_w*TimeHorizon) ];

							K_prefactor{tll_index} = [	S_u{word_id} ;
												eye(n_u*TimeHorizon);
												zeros(tll_card*n_w*TimeHorizon,n_u*TimeHorizon);
												zeros(tll_card*n_x,n_u*TimeHorizon) ] ;

							H{tll_index} = H_PhiI * ( nonK_prefactor{tll_index} ) ;
							
							h_independent_factor{tll_index} = h_PhiI - ...
								H_PhiI * [ J{word_id}*x0 + S_w{word_id}*f_bar{word_id} ;
										zeros(n_u*TimeHorizon,1) ;
										zeros(tll_card*n_w*TimeHorizon,1) ; 
										repmat(x0,tll_card,1) ];

							h_dependent_factor{tll_index} = - ...
								H_PhiI * [ S_u{word_id} ;
										eye(n_u*TimeHorizon) ;
										zeros(tll_card*n_w*TimeHorizon,n_u*TimeHorizon) ; 
										zeros(tll_card*n_x,n_u*TimeHorizon) ];

							h{tll_index} = h_independent_factor{tll_index} + h_dependent_factor{tll_index}*k{knowl_seq_index};

							temp_dummy_w{end+1} = sdpvar(size(H{tll_index},2),1,'full');
							temp_dummy_y{end+1} = sdpvar(size(H{tll_index},1),1,'full');

							% % Create Constraints
							% dummy_var_bound_constraints = dummy_var_bound_constraints + [ -temp_K_bound <= temp_dummy_eta{end} <= temp_K_bound ] + [ temp_dummy_y{end} <= temp_K_bound ];

						end

						% [ reachability_dual_vars{end+1} , incl_constraints ] = cg.get_H_polyt_inclusion_constr( ...
						% 	H{1} , h{1} , ...
						% 	P_target.A*[ zeros(n_x,n_x*TimeHorizon), eye(n_x), zeros(n_x,(n_u+2*n_w)*TimeHorizon+2*n_x) ] , ...
						% 	P_target.b );

						% Create Constraints
						reachability_dual_var_constraints = reachability_dual_var_constraints + [ reachability_dual_vars{end} <= temp_K_bound ];

						% guaranteed_reachability_constraint = guaranteed_reachability_constraint + ...
						% 	iff( matching_behavior{TimeHorizon-1}(knowl_seq_index) == 1 , incl_constraints );
				 	otherwise
				 		error("Unexpected cardinality.")
				end

				% Create Constraints
				if matching_behavior(knowl_seq_index)
					guaranteed_reachability_constraint = guaranteed_reachability_constraint + incl_constraints{end};
				end

			end

		end

		%%%%%%%%%%%%%%%%
		%% Optimize ? %%
		%%%%%%%%%%%%%%%%

		optimization_constraints = 	input_bounds_constraint + ...
									matching_behavior_constraint + ...
									guaranteed_reachability_constraint + causal_detection_constraints %+ ...
									%dummy_var_bound_constraints %+ reachability_dual_var_constraints

		ops = sdpsettings('verbose',1,'debug',1);

		switch solution_approach
			case 'Bilinear (fmincon)'
				optim0 = optimize(optimization_constraints,[],ops);
			case 'Bilinear (Gurobi)'
				ops = sdpsettings(ops,'solver','gurobi');
                % ops.gurobi.Presolve = 0;
                % ops.gurobi.Nonconvex = 2;
				optim0 = optimize(optimization_constraints,[],ops);
			case 'Bilinear (Ipopt)'
				ops = sdpsettings(ops,'solver','ipopt');
				optim0 = optimize(optimization_constraints,[],ops);
            case 'Bilinear (bmibnb)'
                ops = sdpsettings(ops,'solver','bmibnb');
                optim0 = optimize(optimization_constraints,[],ops);
            case 'Bilinear (Baron)'
                ops = sdpsettings(ops,'solver','baron');
                optim0 = optimize(optimization_constraints,[],ops);
            case 'Bilinear (SNOPT)'
                ops = sdpsettings(ops,'solver','snopt');
                optim0 = optimize(optimization_constraints,[],ops);
			case 'Moment Relaxation'
				[ConRelaxed, ~] = momentmodel(optimization_constraints,[],RelaxationOrder)
				optim0 = optimize(ConRelaxed,[],ops);
			otherwise
				error(['Unrecognized solution approach: ' solution_approach ])
		end

		results.Experiment{choice_index}.Runtime = toc(choice_index_start_time);

		results.Experiment{choice_index}.OptimizationOutput = optim0;
		
		selected_final_sequence_indices = matching_behavior;
		results.Experiment{choice_index}.K = {};
		results.Experiment{choice_index}.k = {};
		abbrev_K = {}; abbrev_k = {};
		tempK2 = {};
		for K_index = 1:length(K)
			results.Experiment{choice_index}.K{K_index} = value(K{K_index});
			results.Experiment{choice_index}.k{K_index} = value(k{K_index});
			tempK2{K_index} = value(K{K_index});

			if selected_final_sequence_indices(K_index)
				abbrev_K{end+1} = zeros(n_u*TimeHorizon,n_w*TimeHorizon);
				abbrev_k{end+1} = value(k{K_index});
			end

		end
		results.K2 = tempK2;

		%%%%%%%%%%%%%%%%%%%%%%%
		%% Visualize Results %%
		%%%%%%%%%%%%%%%%%%%%%%%

		%Only visualize if the problem was feasible.

		if optim0.problem == 0 %The problem was feasible.

			temp_valid_L_Sequences = LK_sequences(:,selected_final_sequence_indices);
			lcsas0.X0 = Polyhedron('lb',x0','ub',x0');
			temp_pob_controller = POB_Feedback(lcsas0,abbrev_K,abbrev_k,'PossibleLSequences',temp_valid_L_Sequences);

			results.Experiment{choice_index}.Controller = temp_pob_controller;

			figure;
			hold on;
			plot(P_target,'color','ghostwhite')
			scatter(x0(1),x0(2)) %Plot Initial Condition
			for t = 1:TimeHorizon
				temp_reach = temp_pob_controller.GetReachableSetAt( t , 1 , 'PwT' , LK_w{1} );
				plot(temp_reach,'color','cyan')
			end

			figure;
			hold on;
			plot(P_target,'color','ghostwhite')
			scatter(x0(1),x0(2)) %Plot Initial Condition
			for t = 1:TimeHorizon
				temp_reach = temp_pob_controller.GetReachableSetAt( t , 2 , 'PwT' , LK_w{2} );
				plot(temp_reach,'color','magenta')
			end

		end

	end

	results.OverallRuntime = toc(OverallRuntimeStart);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Determine if Controller Worked %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% x1_T_validation = [ zeros(2,TimeHorizon*2), eye(2) ]*(J1([1:2*(TimeHorizon+1)],[1:2])*x0+S_u1([1:2*(TimeHorizon+1)],[1:TimeHorizon])*results.u([1:TimeHorizon],1));
	% x2_T_validation = [ zeros(2,TimeHorizon*2), eye(2) ]*(J2([1:2*(TimeHorizon+1)],[1:2])*x0+S_u2([1:2*(TimeHorizon+1)],[1:TimeHorizon])*results.u([1:TimeHorizon],2));

	% disp('is x1_T_validation in P_target?')
	% P_target.contains(x1_T_validation)

	% disp('is x2_T_validation in P_target?')
	% P_target.contains(x2_T_validation)

end

function [tf] = knowledge_sequence_contains_prefix( knowl_seq , prefix_seq )
	%Description:
	%	This function attempts to verify whether or not the prefix prefix_seq is a prefix
	% 	of the target knowledge sequence knowl_seq.

	% Constants

	% Algorithm

	for prefix_seq_index = 1:length(prefix_seq)
		if ~(knowl_seq(prefix_seq_index) == prefix_seq(prefix_seq_index))
			tf = false;
			return;
		end
	end

	tf = true;

end

function [ matching_indices, matching_booleans ] = find_all_knowledge_sequences_with_prefix( knowl_seq_matrix , prefix_seq )
	%Description:
	%	This function attempts to verify whether or not the prefix prefix_seq is a prefix
	% 	of the target knowledge sequence knowl_seq.
	%Assumptions:
	%	knowl_seq_matrix has more rows than prefix_seq

	% Constants

	% Algorithm
	matching_booleans = [];

	for knowl_seq_index = 1:size(knowl_seq_matrix,2)
		temp_knowl_seq = knowl_seq_matrix(:,knowl_seq_index);
		matching_booleans(knowl_seq_index) = knowledge_sequence_contains_prefix( temp_knowl_seq , prefix_seq );
	end

	matching_indices = find(matching_booleans);

end

function [H_Phi, h_Phi] = consistent_set_matrices(varargin)
	%consistent_set_matrices.m
	%Description:
	%	[DO NOT TRUST THIS DOCUMENTATION IT IS NOT YET FINISHED]
	%	Finds a polyhedron that describes what pairs of states and input sequences are compatible/feasible from ALL
	%	switching sequences defined by L.
	%
	%	To place the system into clearer focus. We have a Language-Constrained Switched Affine System (LCSAS):
	%
	%	x_{t+1} = A_{q_t} x_t + B_{q_t} u_t + f_{q_t} + w_t
	%	
	%	where:
	%			- q_t is a natural number that describes the current mode at time t
	%			- w_t belongs to the set W_{q_t} which varies with the mode
	%	The consistency set can also be written as"
	%					{ [y]  | }
	%					{ [u]  | }
	%		C(\sigma) = { [w]  | }
	%					{ [v]  | }
	%					{ [x0] | }
	%					{ [x]  | }
	%
	%Inputs:
	%	lcsas 		- An array of Aff_Dyn() objects. Hopefully the dimensions are all appropriately checked so that
	%				  the state is the proper size in all Aff_Dyn(). Inputs as well, etc.
	%	t 			- The time of interest
	%	L 			- The set of words under consideration.
	%				  We would like to find the set of states for which it is possible to reach when under ALL switching
	%				  sequences in this set with the same inputs (albeit with different disturbances).
	%	use_proj 	- Boolean (true or false).
	%				  Used to tell the function to either skip the creation of Consist_set (false) or complete the
	%				  computation of Consist_set (true) which requires projection operations to be called and may be very slow.
	%
	%Example Usage:
	%	[Phi_t_L] = consistent_set(lcsas,t,L)
	%	[Consist_set, full_set] = consistent_set(lcsas,t,L,P_u,P_x0)
	%	[Consist_set, full_set] = consistent_set(lcsas,t,L,P_u,P_x0,'fb_method','state')
	%	[Consist_set, full_set] = consistent_set(lcsas,t,L,P_u,P_x0,'fb_method','state','use_proj',false)
	%
	%Assumptions:
	%	This formulation assumes that the system does not include a disturbed measurements. i.e. We can perfectly observe the state
	%	at each time.

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if nargin < 3
		error('Not enough input arguments.')
	end

	lcsas = varargin{1};
	t = varargin{2};
	L = varargin{3};
	P_u = varargin{4};
	P_x0 = varargin{5};
	
	if ~isa(lcsas,'LCSAS')
		error('Expecting the first input to be a LCSAS object.')
	end

	if ~isa(L,'Language')
		error('Expecting the language input to be a Language object.')
	end

	varargin_idx = 6;
	while varargin_idx <= nargin
		switch varargin{varargin_idx}
			case 'fb_method'
				fb_type = varargin{varargin_idx+1};
				if ~(strcmp(fb_type,'state') || strcmp(fb_type,'output'))
					error(['Invalid feedback type: ' fb_type ])
				end
				varargin_idx = varargin_idx + 2;
			case 'use_proj'
				use_proj = varargin{varargin_idx+1};
				if ~islogical( use_proj )
					error('The flag for ''use_proj'' should be a boolean.')
				end
				varargin_idx = varargin_idx + 2;
			case 'reduce_representation'
				reduce_flag = varargin{varargin_idx+1};
				if ~islogical( reduce_flag )
					error('The flag for ''reduce_flag'' should be a boolean.')
				end
				varargin_idx = varargin_idx + 2;
			otherwise
				error('Unexpected additional input.')
		end
	end

	if (t < 0)
		error(['t must have a value greater than 0.'])
	end

	for word_ind = 1:L.cardinality()
		if t > length(L.words{word_ind})
			error('t should not be larger than any of the words in L.')
		end
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	n_x = size(lcsas.Dyn(1).A,2);
	n_u = size(lcsas.Dyn(1).B,2);
	n_w = size(lcsas.Dyn(1).B_w,2);
	n_y = size(lcsas.Dyn(1).C,1);
	n_v = size(lcsas.Dyn(1).C_v,2);

	if ~exist('fb_type')
		fb_type = 'state';
	end

	if ~exist('use_proj')
		use_proj = true;
	end

	if ~exist('reduce_flag')
		reduce_flag = true;
	end

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	%Create large disturbance set from Cartesian products
	P_wT = 1; P_x0_L = 1;
	for word_idx = 1:length(L.words)
		for symb_idx = 1:t
			P_wT = P_wT * lcsas.Dyn( L.words{word_idx}(symb_idx) ).P_w;
		end
		%Initial Condition Set for Each Word in L
		P_x0_L = P_x0_L * P_x0;
    end
    %Created Disturbance Sets

    P_uT = 1;
    for t_idx = 1:t
        P_uT = P_uT * P_u;
    end

    % Create mpc matrices for each word in the language L
	Hc = {}; Sc = {}; Jc = {}; fc = {}; Cc = {}; Bwc = {}; Cvc = {};
	for word_ind = 1:length(L.words)
		[Hc{word_ind},Sc{word_ind},Cc{word_ind},Jc{word_ind},fc{word_ind},Bwc{word_ind},Cvc{word_ind}] = lcsas.get_mpc_matrices('word',L.words{word_ind}(1:t));
	end

	H_block = []; S_block = []; J_block = []; f_block = [];
	I_blockx = []; I_blockx2 = [];  I_blocky = [];
	C_block = []; Cv_block = [];
	for word_ind = 1:length(L.words)
		H_block(end+[1:size(Hc{word_ind},1)],end+[1:size(Bwc{word_ind},2)]) = Hc{word_ind}*Bwc{word_ind};
		S_block(end+[1:size(Sc{word_ind},1)],[1:size(Sc{word_ind},2)]) = Sc{word_ind};
		J_block(end+[1:size(Jc{word_ind},1)],end+[1:size(Jc{word_ind},2)]) = Jc{word_ind};
		f_block(end+[1:size(Hc{word_ind}*fc{word_ind},1)],1) = Hc{word_ind}*fc{word_ind};

		I_blockx(end+[1:n_x*(t+1)],[1:n_x*(t+1)]) = eye(n_x*(t+1));
		I_blocky(end+[1:n_y*(t+1)],[1:n_y*(t+1)]) = eye(n_y*(t+1));

	end

	% I_block_x0 = []; 
	% block_select_x0 = zeros(L.cardinality()*n_x,size(I_blockx2));
	% for word_ind = 1:L.cardinality()
	% 	block_select_x0(end+[1:n_x],:) = [ zeros(n_x,(n_x*(t+1))*(word_ind-1)) eye(n_x)  ]
	% end

	%% Constructing the Matrices
	if strcmp(fb_type,'state')
	    
		P_eta = P_uT * P_wT * P_x0_L;

		%Create the set of feasible (x,u,w,x0) tuples
		H_Phi = [zeros(size(P_eta.A,1),n_x*(t+1)),P_eta.A];
		H_Phi = [ 	H_Phi ;
					[-I_blockx, S_block, H_block, J_block] ;
					-[-I_blockx, S_block, H_block, J_block] ];

		h_Phi = P_eta.b;
		h_Phi = [ 	h_Phi ;
					-f_block ;
					-(-f_block) ]; %Adding equality constraints

	else

		error('This function is not ready yet.')

		for word_ind = 1:length(L.words)
			C_block(end+[1:n_y*(t+1)],end+[1:n_x*(t+1)]) = [Cc{word_ind} ; zeros(n_y,n_x*t), lcsas.Dyn( L.words{word_ind}(t+1) ).C ];
			Cv_block(end+[1:n_y*(t+1)],end+[1:n_v*(t+1)]) = [Cvc{word_ind},zeros(size(Cvc{word_ind},1),n_v);zeros(n_y,size(Cvc{word_ind},2)), lcsas.Dyn( L.words{word_ind}(t+1) ).C_v ];
			I_blockx2(end+[1:n_x*(t+1)],end+[1:n_x*(t+1)]) = eye(n_x*(t+1));
		end

		%Also introduce the measurement disturbance into the equation
		P_vT = 1; 
		for word_idx = 1:length(L.words)
			for symb_idx = 1:t+1
				P_vT = P_vT * lcsas.Dyn( L.words{word_idx}(symb_idx) ).P_v;
			end
    	end

    	P_eta = P_uT * P_wT * P_vT * P_x0_L;

    	%Create the set of feasible (x,u,w,x0) tuples
    	full_set = Polyhedron(	'A',[zeros(size(P_eta.A,1),n_y*(t+1)),P_eta.A,zeros(size(P_eta.A,1),length(L.words)*n_x*(t+1))],'b',P_eta.b, ...
    							'Ae',[zeros(size(S_block,1),size(I_blocky,2)),S_block, H_block, zeros(size(S_block,1),size(Cv_block,2)), J_block, -I_blockx2; ...
    								  I_blocky, zeros(size(I_blocky,1),size(S_block,2)+size(H_block,2)), -Cv_block , zeros(size(I_blocky,1),size(J_block,2)) , -C_block ], ...
    							'be', [-f_block;zeros(size(I_blocky,1),1)] );

	end

end

function [ possible_LK_realizations , possible_choices ] = get_all_possible_matching_behaviors( lcsas_in , LK )
	%Description:
	%	This function attempts to identify all possible combinations of
	%	knowledge sequences where the knowledge sequences come from LK{end}.

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	num_LK_T_items = length(LK{end});

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	% Identify all possible combinations of matching behaviors that might exist.
	possible_choices = {};
	for temp_length = lcsas_in.L.cardinality():num_LK_T_items
		possible_choices{temp_length} = nchoosek([1:num_LK_T_items],temp_length);
	end

	% possible_choices

	% Translate this choice cell array into a cell array of sequences.

	possible_LK_realizations = {};
	for temp_length = 1:num_LK_T_items
		LK_realizations_at_tl = {};
		for choice_idx = 1:size(possible_choices{temp_length},1)
			temp_choice = possible_choices{temp_length}(choice_idx,:);
			LK_realizations_at_tl{choice_idx} = LK{end}(:,temp_choice);
		end
		possible_LK_realizations{temp_length} = LK_realizations_at_tl;
	end

end

function [ binary_cell_arr ] = convert_closed_loop_LK_to_binary_flags( LK , cl_LK_T )
	%Description:
	%	Identifies if the language sequences from LK are part of the cl_LK.

	%% Algorithm

	binary_cell_arr = {};
	for t = 0:length(LK)-1
		%Create an "all false" array of appropriate dimension
		binary_cell_arr{t+1} = false(length(LK{t+1}),1);

		for seq_index = 1:length(LK{t+1})

			temp_LK_seq = LK{t+1}(:,seq_index);
			seq_i_flag = false;

			for cl_LK_index = 1:size(cl_LK_T,2)
				cl_LK_word_i = cl_LK_T(:,cl_LK_index);
				cl_LK_word_i_prefix = cl_LK_word_i([1:length(temp_LK_seq)]);

				seq_i_flag = seq_i_flag || ( cl_LK_word_i_prefix == temp_LK_seq );

			end

			% Save to the binary_cell_arr
			binary_cell_arr{t+1}(seq_index) = seq_i_flag;

		end
	end

end