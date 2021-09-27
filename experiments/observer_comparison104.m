function [results] = observer_comparison104( varargin )
	%observer_comparison104.m
	%Description:
	%	Attempting to debug why there isn't a differentiating input for the toy2 LCSLS


	disp(' ')
	disp('Beginning observer_comparison104.m')

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

	%% Constants %%

	TimeHorizon = 4;
	[ lcsas , TimeHorizon , Pu , Pw , x0 , Px0 , P_target ] = get_similar_rotation_lcsas('TimeHorizon',TimeHorizon,'eta_u',10);
	X_Target = P_target;


	settings = struct( ...
		'verbosity', 1 , ...
		'subset_search_strategy' , 'AscendingCardinality' , ...
		'UseParallelization' , false , ...
		'DoOptimizationPruningWhere' , 'BeforeSearch' , ...
		'GurobiNodeLimit', 10^7 , ...
		'RemoveBilinearityInReachabilityConstraints', true, ...
		'RemoveBilinearityInInputConstraints', true ...
		);

	settings.DoOptimizationPruningWhere = 'DuringSearch';

	L = lcsas.L;
	TimeHorizon = length(L.words{1});

	[ n_x , n_u , n_y , n_w , n_v ] = lcsas.Dimensions();
	cg = constr_gen(0);

	results.System = lcsas;
	results.Settings = settings;

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	algorithm_start = tic();

	%% Construct All Possible Subsets for this System

	[ unchecked_knowledge_sequences , sequence_construction_history ] = L.create_belief_sequences_of_length(TimeHorizon);
	synthesis_info.UncheckedBeliefSequences = unchecked_knowledge_sequences;

	num_knowl_sequences = size(unchecked_knowledge_sequences,2);

	[ possible_subsets_of_paths , ~ , choices_as_binary_flags ] = lcsas.get_feasible_combinations_of_beliefs( unchecked_knowledge_sequences , ...
																												'verbosity' , settings.verbosity , ... 
																												'SkipOptimizationBasedPruning', strcmp(settings.DoOptimizationPruningWhere,'DuringSearch') );
	synthesis_info.PossibleSubsetsOfPaths = possible_subsets_of_paths;

	[ possible_subsets_of_paths , choices_as_binary_flags ] = organize_subsets_of_paths( unchecked_knowledge_sequences , choices_as_binary_flags , settings.subset_search_strategy );

	synthesis_info.Timing.ConstructingAllSubsets = toc(algorithm_start);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Constructing controller %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	P_search_start = tic();
	synthesis_info.Timing.Experiment = {};

	%for subset_index = 1:length(possible_subsets_of_paths)
	subset_index = 2;

		%Start this iteration (Housekeeping)
		subset_index_start_time = tic;
		active_sequence_flags = choices_as_binary_flags{subset_index};

		%% Make Announcements %%

		if settings.verbosity > 0
			disp(' ')
			disp(['Subset ' num2str(subset_index) '/' num2str(length(possible_subsets_of_paths)) ])
			disp(['- Cardinality of Choice = ' num2str(sum(active_sequence_flags)) ])
			disp(['- Number of Sequences in LK_sequences = ' num2str(size(unchecked_knowledge_sequences,2)) ])
			disp(' ')
		end

		%% Optional Optimization Based Pruning %%

		if strcmp(settings.DoOptimizationPruningWhere,'DuringSearch')
			%Check to see if the proposed paths "cover"
			temp_choice = unchecked_knowledge_sequences(:,active_sequence_flags);

			if ~temp_choice.EvaluateIfAllWordBehaviorsAreCovered( lcsas )
				synthesis_info.About(subset_index).Message = 'Infeasible (System Detected that Subset Did Not Cover)';
				restults.synthesis_info = synthesis_info;
				return;
			end
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

		if settings.RemoveBilinearityInInputConstraints
			strengthen_flag_inputs = 'A_ol';
		else
			strengthen_flag_inputs = 'A_cl';
		end

		for knowl_seq_index = 1:num_knowl_sequences
			knowl_seq = unchecked_knowledge_sequences(:,knowl_seq_index);
			% Create Constraints
			if active_sequence_flags(knowl_seq_index)
			 	ibs_ksi = InternalBehaviorSet(lcsas,knowl_seq, ...
					'OpenLoopOrClosedLoop','Closed',K{knowl_seq_index},k{knowl_seq_index});

			 	[ temp_constraints , input_bound_dual_vars{end+1} ] = ibs_ksi.GetInputBoundConstraints('Use A_cl or A_ol?',strengthen_flag_inputs);

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

			end

	 	end

	 	% Create constraints on which closed loop beliefs ARE NOT POSSIBLE
		dummy_var_bound_constraints = []; temp_dummy_w2 = {}; temp_dummy_y2 = {};
		
		infeasible_sequences = unchecked_knowledge_sequences(:,~active_sequence_flags);
		infeasible_sequences_indices = find(~active_sequence_flags)';
		feasible_sequences = unchecked_knowledge_sequences(:,active_sequence_flags);
		feasible_sequences_indices = find(active_sequence_flags)';

		for infeas_seq_index = 1:size(infeasible_sequences,2)
			infeasible_knowl_seq = infeasible_sequences(:,infeas_seq_index);
			iks_index = infeasible_sequences_indices(infeas_seq_index);

			for feasible_seq_index = 1:size(feasible_sequences,2)
				feasible_knowl_seq = feasible_sequences(:,feasible_seq_index);
				fks_index = feasible_sequences_indices(feasible_seq_index);

				if infeasible_knowl_seq <= feasible_knowl_seq
					%When there is a covering sequence to consider, let's enforce this containment property.
			 		ebs_in = ExternalBehaviorSet( lcsas , infeasible_knowl_seq , ...
			 						'OpenLoopOrClosedLoop','Closed',K{fks_index}, k{fks_index} );

			 		ebs_circum = ExternalBehaviorSet( lcsas , feasible_knowl_seq , ...
			 						'OpenLoopOrClosedLoop','Closed',K{fks_index}, k{fks_index} );

			 		[ temp_dummy_w2{end+1} , temp_constraints ] = CreateContainmentConstraint( ebs_in , ebs_circum );
			 		infeasible_belief_constraints = infeasible_belief_constraints + temp_constraints;
			 	else
			 		ibs_pci = InternalBehaviorSet( lcsas , infeasible_knowl_seq , ...
				  					'OpenLoopOrClosedLoop','Closed',K{fks_index},k{fks_index} );
			 		[ temp_constraints , temp_dummy_y2{end+1} ] = ibs_pci.CreateEmptyConstraint();
			 		infeasible_belief_constraints = infeasible_belief_constraints + temp_constraints;
			 	end
			 		
			end

	 	end


	 	% Create Causal Detection Constraints %%
		causal_detection_constraints = cg.get_belief_prefix_gain_constraints( lcsas , K , k , unchecked_knowledge_sequences );
		causal_gain_constraints = lcsas.create_lower_diagonal_constraint_on_gains( K , 'Disturbance' );

		% Create Robust Reachability Constraints

		guaranteed_reachability_constraint = [];
		reachability_dual_vars = {};

		if settings.RemoveBilinearityInReachabilityConstraints
			strengthen_flag_reachability = 'A_ol';
		else
			strengthen_flag_reachability = 'A_cl';
		end
        
		for knowl_seq_index = 1:num_knowl_sequences
			
			knowl_seq = unchecked_knowledge_sequences(:,knowl_seq_index);
			
			% Create Constraints
			if active_sequence_flags(knowl_seq_index)

				ibs_ksi = InternalBehaviorSet( lcsas , knowl_seq , ...
												'OpenLoopOrClosedLoop','Closed', K{knowl_seq_index}, k{knowl_seq_index});
				
				[ temp_constraints , reachability_dual_vars{end+1} ] = ibs_ksi.GetReachabilityConstraints( X_Target , 'Use A_cl or A_ol?' , strengthen_flag_reachability );

				guaranteed_reachability_constraint = guaranteed_reachability_constraint + temp_constraints;
			end

		end

		% Optimize!

		optimization_constraints = 	feasible_belief_constraints + infeasible_belief_constraints + ...
									input_bounds_constraints + ...
									guaranteed_reachability_constraint + causal_detection_constraints + causal_gain_constraints

		ops = sdpsettings(	'verbose',settings.verbosity,'debug',1);
		ops.gurobi.NodeLimit = settings.GurobiNodeLimit;

		ops = sdpsettings(ops,'solver','gurobi');
		optim0 = optimize(optimization_constraints,[],ops);

		synthesis_info.Timing.Experiment{subset_index}.Optimization.yalmiptime = optim0.yalmiptime;
		synthesis_info.Timing.Experiment{subset_index}.Optimization.solvertime = optim0.solvertime; 
		synthesis_info.Timing.Experiment{subset_index}.ForExperiment = toc(subset_index_start_time);

		controller = [];
		if optim0.problem == 0
			%% Create Controller
			
			synthesis_info.K = {};
			synthesis_info.k = {};

			for K_index = 1:length(K)
				synthesis_info.K{K_index} = value(K{K_index});
				synthesis_info.k{K_index} = value(k{K_index});
			end

			for Lambda_index = 1:length(reachability_dual_vars)
				LambdaReach{Lambda_index} = value(reachability_dual_vars{Lambda_index});
			end
			synthesis_info.LambdaReach = LambdaReach

			controller = ConsistentBeliefsController( lcsas , unchecked_knowledge_sequences(:,active_sequence_flags) , ...
														{ synthesis_info.K{active_sequence_flags} } , ...
														{ synthesis_info.k{active_sequence_flags} } );

			[ tf , norm_matrix_diff , vector_diff ] = check_reachability_condition( controller , LambdaReach , X_Target , ...
																					settings.RemoveBilinearityInReachabilityConstraints );
			synthesis_info.ReachabilityConstraintData.Satisfied = tf;
			synthesis_info.ReachabilityConstraintData.NormMatrixDiff = norm_matrix_diff;
			synthesis_info.ReachabilityConstraintData.VectorDiff = vector_diff;

			synthesis_info.About(subset_index).problem = optim0.problem;
			synthesis_info.About(subset_index).Message = 'Solved!';

			results.synthesis_info = synthesis_info;

			return;
		else
			synthesis_info.About(subset_index).problem = optim0.problem;
			synthesis_info.About(subset_index).Message = 'Not Solved (Optimization''s ''problem'' field was not satisfactory.)';
		end


	% end

	synthesis_info.Timing.SearchForFeasibleP = toc(P_search_start);
	synthesis_info.Timing.Overall = toc(algorithm_start);

	results.synthesis_info = synthesis_info;

end