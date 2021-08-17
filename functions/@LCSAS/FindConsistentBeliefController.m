function [ controller , synthesis_info ] = FindConsistentBeliefController( varargin )
	%Description:
	%	Finds a ConsistentBeliefController according to the method defined in our Journal submission.
	%
	%Usage:
	%	[ cbc_out , synthesis_info ] = system.FindConsistentBeliefController( X_Target )

	%% Input Processing %%

	[ lcsas , X_Target , settings ] = ip_FindConsistentBeliefController( varargin{:} );

	synthesis_info.System = lcsas;
	synthesis_info.TargetSet = X_Target;

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	L = lcsas.L;
	TimeHorizon = length(L.words{1});

	[ n_x , n_u , n_y , n_w , n_v ] = lcsas.Dimensions();
	cg = constr_gen(0);

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	algorithm_start = tic();

	%% Construct All Possible Subsets for this System

	[ unchecked_knowledge_sequences , sequence_construction_history ] = L.create_belief_sequences_of_length(TimeHorizon);
	synthesis_info.UncheckedBeliefSequences = unchecked_knowledge_sequences;

	num_knowl_sequences = size(unchecked_knowledge_sequences,2);

	[ possible_subsets_of_paths , ~ , choices_as_binary_flags ] = lcsas.get_feasible_combinations_of_beliefs( unchecked_knowledge_sequences , 'verbosity' , settings.verbosity );
	synthesis_info.PossibleSubsetsOfPaths = possible_subsets_of_paths;

	[ possible_subsets_of_paths , choices_as_binary_flags ] = organize_subsets_of_paths( unchecked_knowledge_sequences , choices_as_binary_flags , settings.subset_search_strategy );

	synthesis_info.Timing.ConstructingAllSubsets = toc(algorithm_start);

	%% Constructing controller

	P_search_start = tic();
	synthesis_info.Timing.Experiment = {};

	for subset_index = 1:length(possible_subsets_of_paths)

		%Start this iteration
		subset_index_start_time = tic;
		active_sequence_flags = choices_as_binary_flags{subset_index};

		if settings.verbosity > 0
			disp(' ')
			disp(['Subset ' num2str(subset_index) '/' num2str(length(possible_subsets_of_paths)) ])
			disp(['- Cardinality of Choice = ' num2str(sum(active_sequence_flags)) ])
			disp(['- Number of Sequences in LK_sequences = ' num2str(size(unchecked_knowledge_sequences,2)) ])
			disp(' ')
		end


		% Create Optimization Variables
		K = {}; k = {};
		for knowl_seq_index = 1:num_knowl_sequences
			%K{knowl_seq_index} = zeros(TimeHorizon*n_u,TimeHorizon*n_w);
			K{knowl_seq_index} = sdpvar(TimeHorizon*n_u,(TimeHorizon)*n_w,'full');
			k{knowl_seq_index} = sdpvar(TimeHorizon*n_u,1,'full');
		end

		% Input Bound Constraints
		input_bounds_constraints = [];
		input_bound_dual_vars = {};

		for knowl_seq_index = 1:num_knowl_sequences
			knowl_seq = unchecked_knowledge_sequences(:,knowl_seq_index);
			% Create Constraints
			if active_sequence_flags(knowl_seq_index)
			 	ibs_ksi = InternalBehaviorSet(lcsas,knowl_seq, ...
					'OpenLoopOrClosedLoop','Closed',K{knowl_seq_index},k{knowl_seq_index});

			 	[ temp_constraints , input_bound_dual_vars{end+1} ] = ibs_ksi.GetInputBoundConstraints('Relaxation',true);

			 	input_bounds_constraints = input_bounds_constraints + temp_constraints;
			 end 
		end

		% Matching Behavior Constraints

		feasible_belief_constraints = [];
		infeasible_belief_constraints = [];

		active_sequence_flags_constraint = [];
		dummy_var_bound_constraints = []; temp_dummy_w = {}; temp_dummy_y = {};
		for knowl_seq_index = 1:num_knowl_sequences
			
			knowl_seq = unchecked_knowledge_sequences(:,knowl_seq_index);
			temp_last_lang = knowl_seq(end);
			tll_card = temp_last_lang.cardinality();

			% Create Constraints
			if active_sequence_flags(knowl_seq_index)

				ibs_ksi = InternalBehaviorSet(lcsas,knowl_seq, ...
					'OpenLoopOrClosedLoop','Closed',K{knowl_seq_index},k{knowl_seq_index});

			 	[ temp_constraints, temp_dummy_w{end+1} ] = ibs_ksi.CreateNonemptyConstraint();

			 	feasible_belief_constraints = feasible_belief_constraints + temp_constraints;

			 	% Construct P^C( knowl_seq )
			 	P_C_indices = ~active_sequence_flags;
			 	for pci_index = 1:length(P_C_indices)
			 		%If this sequence is in the complement
			 		if P_C_indices(pci_index)
			 			%Grab path
			 			temp_seq = unchecked_knowledge_sequences(:,pci_index);
			 			%If this sequence IS NOT a subset of knowl_seq, then enforce empty constraint of closed loop.
			 			if ~temp_seq.contains(knowl_seq)
			 				ibs_pci = InternalBehaviorSet(lcsas,temp_seq, ...
			 					'OpenLoopOrClosedLoop','Closed',K{knowl_seq_index},k{knowl_seq_index});

			 				[ temp_constraints , temp_dummy_y{end+1} ] = ibs_pci.CreateEmptyConstraint();
			 				infeasible_belief_constraints = infeasible_belief_constraints + temp_constraints;
			 			end

			 		end
			 	end

			end

	 	end

	 	% Create Causal Detection Constraints %%
		causal_detection_constraints = cg.get_belief_prefix_gain_constraints( lcsas , K , k , unchecked_knowledge_sequences );
		causal_gain_constraints = lcsas.create_lower_diagonal_constraint_on_gains( K , 'Disturbance' );

		% Create Robust Reachability Constraints

		guaranteed_reachability_constraint = [];
		reachability_dual_vars = {};
        
		for knowl_seq_index = 1:num_knowl_sequences
			
			knowl_seq = unchecked_knowledge_sequences(:,knowl_seq_index);
			
			% Create Constraints
			if active_sequence_flags(knowl_seq_index)

				ibs_ksi = InternalBehaviorSet( lcsas , knowl_seq , ...
												'OpenLoopOrClosedLoop','Closed', K{knowl_seq_index}, k{knowl_seq_index});
				
				[ temp_constraints , reachability_dual_vars{end+1} ] = ibs_ksi.GetReachabilityConstraints( X_Target , 'Relaxation' , true );

				guaranteed_reachability_constraint = guaranteed_reachability_constraint + temp_constraints;
			end

		end

		% Optimize!

		optimization_constraints = 	feasible_belief_constraints + infeasible_belief_constraints + ...
									input_bounds_constraints + ...
									guaranteed_reachability_constraint + causal_detection_constraints + causal_gain_constraints

		ops = sdpsettings(	'verbose',settings.verbosity,'debug',1);
		ops.gurobi.NodeLimit = 10^5;

		ops = sdpsettings(ops,'solver','gurobi');
		optim0 = optimize(optimization_constraints,[],ops);

		synthesis_info.Timing.Experiment{subset_index}.Optimization.yalmiptime = optim0.yalmiptime;
		synthesis_info.Timing.Experiment{subset_index}.Optimization.solvertime = optim0.solvertime; 
		synthesis_info.Timing.Experiment{subset_index}.ForExperiment = toc(subset_index_start_time);

		if optim0.problem == 0
			break;
		end


	end

	synthesis_info.Timing.SearchForFeasibleP = toc(P_search_start);
	synthesis_info.Timing.Overall = toc(algorithm_start);

	controller = [];

end

function [ lcsas , X_Target , settings ] = ip_FindConsistentBeliefController( varargin )
	%Description:
	%	Processing the inputs that are given to the FindConsistentBeliefController

	%% Constants %%

	% Default Settings
	settings = struct( ...
		'verbosity', 1 , ...
		'subset_search_strategy' , 'AscendingCardinality' ...
		);

	%% Algorithm %%

	lcsas = varargin{1};
	X_Target = varargin{2};

	if nargin > 2
		argument_index = 3;
		while argument_index <= nargin
			switch varargin{argument_index}
				case 'SearchStrategy'
					settings.subset_search_strategy = varargin{argument_index+1};
					argument_index = argument_index + 2;
				otherwise
					error(['Unexpected input to FindConsistentBeliefController: ' varargin{argument_index} ]);
			end
		end
	end

end