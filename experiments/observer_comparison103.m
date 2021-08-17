function [results] = observer_comparison103( varargin )
	%observer_comparison103.m
	%Description:
	%	Attempting to implement a solution to subproblem 3 of the LCSLS wiki pages
	%	using an enumeration based version of the approach given in observer_comparison86.m.
	%	This method also restricts itself by not using linear gains.

	disp(' ')
	disp('Beginning observer_comparison102.m')
	disp('Directly running the enumeration approach for similar rotation system')
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

% 			lcsas0.L = Language(2*ones(1,TimeHorizon));
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

	[ possibleSubsetsOfPaths , possible_choices , choices_as_binary_flags ] = lcsas0.get_feasible_combinations_of_beliefs( LK_sequences , 'verbosity' , 1 );

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

	%choice_index = length(choices_as_binary_flags);
	for choice_index = 1:length(choices_as_binary_flags)

		disp(' ')
		disp(['Choice ' num2str(choice_index) '/' num2str(length(choices_as_binary_flags)) ])
		disp(['- Cardinality of Choice = ' num2str(sum(choices_as_binary_flags{choice_index})) ])
		disp(['- Number of Sequences in LK_sequences = ' num2str(size(LK_sequences,2)) ])
		disp(' ')

		%Start timer
		choice_index_start_time = tic;

		% Variables

		% Create Optimization Variables
		K = {}; k = {};
		for knowl_seq_index = 1:num_knowl_sequences
			%K{knowl_seq_index} = zeros(TimeHorizon*n_u,TimeHorizon*n_w);
			K{knowl_seq_index} = sdpvar(TimeHorizon*n_u,(TimeHorizon)*n_w,'full');
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
		input_bounds_constraints = [];
		input_bound_dual_vars = {};

		for knowl_seq_index = 1:num_knowl_sequences
			knowl_seq = LK_sequences(:,knowl_seq_index);
			% Create Constraints
			if matching_behavior(knowl_seq_index)
			 	ibs_ksi = InternalBehaviorSet(lcsas0,knowl_seq, ...
					'OpenLoopOrClosedLoop','Closed',K{knowl_seq_index},k{knowl_seq_index});

			 	[ temp_constraints , input_bound_dual_vars{end+1} ] = ibs_ksi.GetInputBoundConstraints('Relaxation',true);

			 	input_bounds_constraints = input_bounds_constraints + temp_constraints;
			 end 
		end

		% input_bounds_constraint = cg.get_input_bound_constraint_on( lcsas0 , K , k , 'fb_type' , 'state-disturbance' ); 

		%% Matching Behavior Constraints

		feasible_belief_constraints = [];
		infeasible_belief_constraints = [];

		matching_behavior_constraint = [];
		dummy_var_bound_constraints = []; temp_dummy_w = {}; temp_dummy_y = {};
		for knowl_seq_index = 1:num_knowl_sequences
			
			knowl_seq = LK_sequences(:,knowl_seq_index);
			temp_last_lang = knowl_seq(end);
			tll_card = temp_last_lang.cardinality();

			% Create Constraints
			if matching_behavior(knowl_seq_index)

				ibs_ksi = InternalBehaviorSet(lcsas0,knowl_seq, ...
					'OpenLoopOrClosedLoop','Closed',K{knowl_seq_index},k{knowl_seq_index});

			 	[ temp_constraints, temp_dummy_w{end+1} ] = ibs_ksi.CreateNonemptyConstraint();

			 	feasible_belief_constraints = feasible_belief_constraints + temp_constraints;

			 	% Construct P^C( knowl_seq )
			 	P_C_indices = ~matching_behavior;
			 	for pci_index = 1:length(P_C_indices)
			 		%If this sequence is in the complement
			 		if P_C_indices(pci_index)
			 			%Grab path
			 			temp_seq = LK_sequences(:,pci_index);
			 			%If this sequence IS NOT a subset of knowl_seq, then enforce empty constraint of closed loop.
			 			if ~temp_seq.contains(knowl_seq)
			 				ibs_pci = InternalBehaviorSet(lcsas0,temp_seq, ...
			 					'OpenLoopOrClosedLoop','Closed',K{knowl_seq_index},k{knowl_seq_index});

			 				[ temp_constraints , temp_dummy_y{end+1} ] = ibs_pci.CreateEmptyConstraint();
			 				infeasible_belief_constraints = infeasible_belief_constraints + temp_constraints;
			 			end

			 		end
			 	end

			end

		end

		%% Create Causal Detection Constraints %%
		causal_detection_constraints = cg.get_belief_prefix_gain_constraints( lcsas0 , K , k , LK_sequences );
		causal_gain_constraints = lcsas0.create_lower_diagonal_constraint_on_gains( K , 'Disturbance' );

		%% Create Robust Reachability Constraints

		guaranteed_reachability_constraint = [];
		reachability_dual_vars = {};
        
		for knowl_seq_index = 1:num_knowl_sequences
			
			knowl_seq = LK_sequences(:,knowl_seq_index);
			
			% Create Constraints
			if matching_behavior(knowl_seq_index)

				ibs_ksi = InternalBehaviorSet( lcsas0 , knowl_seq , ...
												'OpenLoopOrClosedLoop','Closed', K{knowl_seq_index}, k{knowl_seq_index});
				
				[ temp_constraints , reachability_dual_vars{end+1} ] = ibs_ksi.GetReachabilityConstraints( P_target , 'Relaxation' , true );

				guaranteed_reachability_constraint = guaranteed_reachability_constraint + temp_constraints;
			end

		end

		%%%%%%%%%%%%%%%%
		%% Optimize ? %%
		%%%%%%%%%%%%%%%%

		optimization_constraints = 	feasible_belief_constraints + infeasible_belief_constraints + ...
									input_bounds_constraints + matching_behavior_constraint + ...
									guaranteed_reachability_constraint + causal_detection_constraints + causal_gain_constraints %+ ...
		
		%optimization_constraints = 	feasible_belief_constraints ; % SUCCESSFUL
		%optimization_constraints = 	feasible_belief_constraints + infeasible_belief_constraints; % UNSUCCESSFUL
		%optimization_constraints = 	feasible_belief_constraints + input_bounds_constraints; %SUCCESSFUL
		%optimization_constraints = 	feasible_belief_constraints + input_bounds_constraints + matching_behavior_constraint; %SUCCESSFUL
		%optimization_constraints = 	feasible_belief_constraints + input_bounds_constraints + matching_behavior_constraint + guaranteed_reachability_constraint; %UNSUCCESSFUL

		ops = sdpsettings(	'verbose',1,'debug',1);
		ops.gurobi.NodeLimit = 10^5;

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

		for Lambda_index = 1:length(reachability_dual_vars)
			LambdaReach{Lambda_index} = value(reachability_dual_vars{Lambda_index});
		end
		results.Experiment{choice_index}.LambdaReach = LambdaReach

		%%%%%%%%%%%%%%%%%%%%%%%
		%% Visualize Results %%
		%%%%%%%%%%%%%%%%%%%%%%%

		%Only visualize if the problem was feasible.

		results.Experiment{choice_index}.SimulationData = [];

		if optim0.problem == 0

			cbc0 = ConsistentBeliefsController( lcsas0 , LK_sequences(:,matching_behavior) , ...
												results.Experiment{choice_index}.K, ...
												results.Experiment{choice_index}.k );

			[ tf , norm_matrix_diff , vector_diff ] = check_reachability_condition( cbc0 , LambdaReach , P_target , true )
			results.Experiment{choice_index}.ReachabilityConstraintData.Satisfied = tf;
			results.Experiment{choice_index}.ReachabilityConstraintData.NormMatrixDiff = norm_matrix_diff;
			results.Experiment{choice_index}.ReachabilityConstraintData.VectorDiff = vector_diff;

			figure;

			hold on;
			plot(lcsas0.X0)
			plot(P_target)

			for simulation_index = 1:10
				[ x_0_t, u_0_tm1 , y_0_t , sig ] = cbc0.simulate_1run();
				% for t = 0:TimeHorizon
				% 	x_t = x_0_t(:,t+1);
				% 	scatter(x_t(1),x_t(2))
				% end
				plot(x_0_t(1,:),x_0_t(2,:))
				cbc0.clear_histories()

				%Save data
				results.Experiment{choice_index}.SimulationData = [ results.Experiment{choice_index}.SimulationData ; struct('x_0_t',x_0_t,'u_0_tm1',u_0_tm1) ];

			end

			axis([-3,4,-4,4])


			return;
		end



		% if optim0.problem == 0 %The problem was feasible.

		% 	temp_valid_L_Sequences = LK_sequences(:,selected_final_sequence_indices);
		% 	lcsas0.X0 = Polyhedron('lb',x0','ub',x0');
		% 	temp_pob_controller = POB_Feedback(lcsas0,abbrev_K,abbrev_k,'PossibleLSequences',temp_valid_L_Sequences);

		% 	results.Experiment{choice_index}.Controller = temp_pob_controller;

		% 	figure;
		% 	hold on;
		% 	plot(P_target,'color','ghostwhite')
		% 	scatter(x0(1),x0(2)) %Plot Initial Condition
		% 	for t = 1:TimeHorizon
		% 		temp_reach = temp_pob_controller.GetReachableSetAt( t , 1 , 'PwT' , LK_w{1} );
		% 		plot(temp_reach,'color','cyan')
		% 	end

		% 	figure;
		% 	hold on;
		% 	plot(P_target,'color','ghostwhite')
		% 	scatter(x0(1),x0(2)) %Plot Initial Condition
		% 	for t = 1:TimeHorizon
		% 		temp_reach = temp_pob_controller.GetReachableSetAt( t , 2 , 'PwT' , LK_w{2} );
		% 		plot(temp_reach,'color','magenta')
		% 	end

		% end

		choice_index = choice_index - 1;

	end

	results.OverallRuntime = toc(OverallRuntimeStart);

end