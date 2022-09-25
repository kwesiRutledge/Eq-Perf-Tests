function [ controller , synthesis_info ] = FindAdaptiveController_AlwaysHidenInfo( varargin )
	%Description:
	%	Finds a ConsistentBeliefController according to the method defined in our Journal submission.
	%
	%Inputs:
	%	- X_Target: Polytope representing the set of states where the 
	%	- 'SearchStrategy': Optional input argument.
	%						Defines how the algorithm will search through the set of all possible KnowledgeSequence objects (i.e. 
	%						P-tilde objects in the math). Can be in 'DescendingCardinality', 'AscendingCardinality', ...
	%						By default the search strategy is 'AscendingCardinality'
	%	-
	%
	%Usage:
	%	[ cbc_out , synthesis_info ] = system.FindConsistentBeliefController( X_Target )
	%	[ cbc_out , synthesis_info ] = system.FindConsistentBeliefController( X_Target , 'SearchStrategy' , 'DescendingCardinality' )

	%% Input Processing %%

	[ lcsas , X_Target , settings ] = ip_FindAdaptiveController_AlwaysHidenInfo( varargin{:} );

	synthesis_info.System = lcsas;
	synthesis_info.TargetSet = X_Target;
	synthesis_info.SynthesisSettings = settings;

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

	controller = [];
	synthesis_info = [];

	algorithm_start = tic();

	%% Construct All Possible Subsets for this System

	[ M_unchecked , sequence_construction_history ] = L.create_belief_sequences_of_length(TimeHorizon);
    M_target = [L];
    for row_index = 2:TimeHorizon
        M_target(row_index,1) = L;
    end
	synthesis_info.UncheckedBeliefSequences = M_unchecked;
    synthesis_info.M_target = M_target;

	synthesis_info.Timing.ConstructMStar = toc(algorithm_start);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Constructing controller %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	P_search_start = tic();
	synthesis_info.Timing.Experiment = {};

	%Start this iteration (Housekeeping)
	subset_index_start_time = tic;

	num_unchecked = size(M_unchecked,2);
    col_index_of_M_target = num_unchecked;

	%% Make Announcements %%

	if settings.verbosity > 0
		disp(' ')
		disp(['M_target for Hiding Info'])
		disp(['- Cardinality of M_target = ' num2str(size(M_target,2)) ])
		%disp(['- Number of Sequences in LK_sequences = ' num2str(size(MStar,2)) ])
		disp(' ')
	end

	% Create Optimization Variables
	K = {}; k = {};
	for m_index = 1:num_unchecked
		%K{knowl_seq_index} = zeros(TimeHorizon*n_u,TimeHorizon*n_w);
		K{m_index} = sdpvar(TimeHorizon*n_u,(TimeHorizon)*n_w,'full');
		k{m_index} = sdpvar(TimeHorizon*n_u,1,'full');
	end

	% Input Bound Constraints
	input_bounds_constraints = [];
	input_bound_dual_vars = {};

	if settings.RemoveBilinearityInInputConstraints
		strengthen_flag_inputs = 'A_ol';
	else
		strengthen_flag_inputs = 'A_cl';
	end

	for m_index = 1:num_unchecked

        %Skip all gains except the one for M_target
        % if m_index ~= num_unchecked
        %     continue;
        % end

		m = M_unchecked(:,m_index);
		% Create Constraints
	 	ibs_m = InternalBehaviorSet(lcsas,m, ...
			'OpenLoopOrClosedLoop','Closed',K{m_index},k{m_index});

	 	[ temp_constraints , input_bound_dual_vars{end+1} ] = ibs_m.GetInputBoundConstraints('Use A_cl or A_ol?',strengthen_flag_inputs);

	 	input_bounds_constraints = input_bounds_constraints + temp_constraints;
	end

	% Matching Behavior Constraints

	feasible_belief_constraints = [];
	infeasible_belief_constraints = [];

	active_sequence_flags_constraint = [];
	dummy_var_bound_constraints = []; temp_dummy_w = {}; temp_dummy_y = {};

 	% Create constraints on which closed loop beliefs ARE NOT POSSIBLE
	dummy_var_bound_constraints = []; temp_dummy_w2 = {}; temp_dummy_y2 = {};

	for m_c_index = 1:num_unchecked
		m_c = M_unchecked(:,m_c_index);
% 		iks_index = infeasible_sequences_indices(infeas_seq_index);
        
        R_mc_prime = ExternalBehaviorSet( lcsas , m_c );

        if R_mc_prime.ToPolyhedron().isEmptySet()
            continue;
        end

		for m_index = num_unchecked
			m = M_unchecked(:,m_index);
			
            % Skip any combination where m == m_c
            if m_c == m
                continue
            end

			if m_c <= m
				%When there is a covering sequence to consider, let's enforce this containment property.
		 		R_mc = ExternalBehaviorSet( lcsas , m_c , ...
		 						'OpenLoopOrClosedLoop','Closed',K{m_index}, k{m_index} );

		 		% ebs_circum = ExternalBehaviorSet( lcsas , feasible_knowl_seq , ...
		 		% 				'OpenLoopOrClosedLoop','Closed', K{fks_index}, k{fks_index} );
		 		R_m = ExternalBehaviorSet( lcsas , m );


		 		temp_constraints = [];
		 		if settings.LinearizeBilinearContainment
		 			[ temp_dummy_w2{end+1} , temp_constraints ] = CreateLinearizedContainmentConstraint( R_mc , R_m );
		 		else
		 			[ temp_dummy_w2{end+1} , temp_constraints ] = CreateContainmentConstraint( R_mc , R_m );
		 		end


		 		infeasible_belief_constraints = infeasible_belief_constraints + temp_constraints;
		 	else
		 		ibs_mc = InternalBehaviorSet( lcsas , m_c , ...
			  					'OpenLoopOrClosedLoop','Closed',K{m_index},k{m_index} );
		 		[ temp_constraints , temp_dummy_y2{end+1} ] = ibs_mc.CreateEmptyConstraint();
		 		infeasible_belief_constraints = infeasible_belief_constraints + temp_constraints;
		 	end

		 		
		end

 	end


 	% Create Causal Detection Constraints %%
	causal_detection_constraints = cg.get_belief_prefix_gain_constraints( lcsas , K , k , M_unchecked );
	causal_gain_constraints = lcsas.create_lower_diagonal_constraint_on_gains( K , 'Disturbance' );

	% Create Robust Reachability Constraints

	guaranteed_reachability_constraint = [];
	reachability_dual_vars = {};

	if settings.RemoveBilinearityInReachabilityConstraints
		strengthen_flag_reachability = 'A_ol';
	else
		strengthen_flag_reachability = 'A_cl';
	end
    
	for m_index = 1:num_unchecked

        %Skip all m that do not belong to M_target
        if m_index ~= num_unchecked
            continue;
        end
		
		m = M_unchecked(:,m_index);
		
		% Create Constraints
		ibs_ksi = InternalBehaviorSet( lcsas , m , ...
										'OpenLoopOrClosedLoop','Closed', K{m_index}, k{m_index});
		X_Target
		[ temp_constraints , reachability_dual_vars{end+1} ] = ibs_ksi.GetReachabilityConstraints( X_Target , 'Use A_cl or A_ol?' , strengthen_flag_reachability );

		guaranteed_reachability_constraint = guaranteed_reachability_constraint + temp_constraints;

	end

	synthesis_info.Timing.ConstraintPrepTime = toc(subset_index_start_time);

	% Optimize!

	optimization_constraints = 	feasible_belief_constraints + infeasible_belief_constraints + ...
								input_bounds_constraints + ...
								guaranteed_reachability_constraint + causal_detection_constraints + causal_gain_constraints;

	ops = sdpsettings(	'verbose',settings.verbosity,'debug',1);
	ops.gurobi.NodeLimit = settings.GurobiNodeLimit;

	ops = sdpsettings(ops,'solver','gurobi');
	optim0 = optimize(optimization_constraints,[],ops);

	synthesis_info.Timing.Optimization.yalmiptime = optim0.yalmiptime;
	synthesis_info.Timing.Optimization.solvertime = optim0.solvertime; 
	synthesis_info.Timing.ForMStar = toc(subset_index_start_time);

	controller = [];

	if optim0.problem == 0
		%% Create Controller
		
		synthesis_info.K = {};
		synthesis_info.k = {};
% 
% 		for K_index = 1:length(K)
			synthesis_info.K{1} = value(K{end});
			synthesis_info.k{1} = value(k{end});
% 		end

% 		for Lambda_index = 1:length(reachability_dual_vars)
% 			LambdaReach{Lambda_index} = value(reachability_dual_vars{Lambda_index});
% 		end
% 		synthesis_info.LambdaReach = LambdaReach;

		controller = ConsistentBeliefsController( lcsas , M_target , synthesis_info.K , synthesis_info.k );
% 
% 		[ tf , norm_matrix_diff , vector_diff ] = check_reachability_condition( controller , LambdaReach , X_Target , ...
% 																				settings.RemoveBilinearityInReachabilityConstraints );
% 		synthesis_info.ReachabilityConstraintData.Satisfied = tf;
% 		synthesis_info.ReachabilityConstraintData.NormMatrixDiff = norm_matrix_diff;
% 		synthesis_info.ReachabilityConstraintData.VectorDiff = vector_diff;

		synthesis_info.Optimization.problem = optim0.problem;
		synthesis_info.Message = 'Solved!';

	else
		synthesis_info.Optimization.problem = optim0.problem;
		synthesis_info.Message = 'Not Solved (Optimization''s ''problem'' field was not satisfactory.)';
	end

	synthesis_info.Timing.SearchForFeasibleP = toc(P_search_start);
	synthesis_info.Timing.Overall = toc(algorithm_start);

end

function [ lcsas , X_Target , settings ] = ip_FindAdaptiveController_AlwaysHidenInfo( varargin )
	%Description:
	%	Processing the inputs that are given to the FindConsistentBeliefController

	%% Constants %%

	% Default Settings
	settings = struct( ...
		'verbosity', 1 , ...
		'subset_search_strategy' , 'AscendingCardinality' , ...
		'UseParallelization' , false , ...
		'DoOptimizationPruningWhere' , 'BeforeSearch' , ...
		'GurobiNodeLimit', 10^5 , ...
		'RemoveBilinearityInReachabilityConstraints', false, ...
		'RemoveBilinearityInInputConstraints', true , ...
		'LinearizeBilinearContainment', false , ...
		'AlgorithmVersion',2 ... 								%Version 1 of the algorithm considers a nonemptiness constraint on each knowledge sequence
		...														%Version 2 does not consider the nonemptiness constraint.
		);

	PruneBeliefsWhere_options = {'BeforeSearch','DuringSearch','Nowhere'};

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
				case 'UseParallelization'
					settings.UseParallelization = varargin{argument_index+1};
					argument_index = argument_index + 2;
				case 'DoOptimizationPruningWhere'
					settings.DoOptimizationPruningWhere = varargin{argument_index+1};
					if ~any(strcmp(settings.DoOptimizationPruningWhere,PruneBeliefsWhere_options))
						error(['The input ' varargin{argument_index+1} ' for DoOptimizationPruningWhere is not allowed!'])
					end
					argument_index = argument_index + 2;
				case 'GurobiNodeLimit'
					settings.GurobiNodeLimit = varargin{argument_index+1};
					argument_index = argument_index + 2;
				case 'RemoveBilinearityInReachabilityConstraints'
					settings.RemoveBilinearityInReachabilityConstraints = varargin{argument_index+1};
					argument_index = argument_index + 2;
				case 'RemoveBilinearityInInputConstraints'
					settings.RemoveBilinearityInInputConstraints = varargin{argument_index+1};
					argument_index = argument_index + 2;
				case 'LinearizeBilinearContainment'
					settings.LinearizeBilinearContainment = varargin{argument_index+1};
					argument_index = argument_index + 2;
				case 'AlgorithmVersion'
					settings.AlgorithmVersion = varargin{argument_index+1};
					argument_index = argument_index + 2;
				otherwise
					error(['Unexpected input to FindConsistentBeliefController: ' varargin{argument_index} ]);
			end
		end
	end

end

function [ controller , synthesis_info ] = FindConsistentBeliefController_Parallelized( lcsas , X_Target , settings )

	%% Constants

	L = lcsas.L;
	TimeHorizon = length(L.words{1});

	[ n_x , n_u , n_y , n_w , n_v ] = lcsas.Dimensions();
	cg = constr_gen(0);

	p = gcp(); %Get Parameters of the current parallel pool.
	p.NumWorkers;

	%% Algorithm

	algorithm_start = tic();

	%% Construct All Possible Subsets for this System

	[ MStar , sequence_construction_history ] = L.create_belief_sequences_of_length(TimeHorizon);
	synthesis_info.UncheckedBeliefSequences = MStar;

	num_knowl_sequences = size(MStar,2);

	[ possible_subsets_of_paths , ~ , choices_as_binary_flags ] = lcsas.get_feasible_combinations_of_beliefs( MStar , 'verbosity' , settings.verbosity );
	synthesis_info.PossibleSubsetsOfPaths = possible_subsets_of_paths;

	[ possible_subsets_of_paths , choices_as_binary_flags ] = organize_subsets_of_paths( MStar , choices_as_binary_flags , settings.subset_search_strategy );

	synthesis_info.Timing.ConstructingAllSubsets = toc(algorithm_start);

	%% Constructing controller

	P_search_start = tic();
	synthesis_info.Timing.Experiment = {};

	F = {};
	for subset_index = 1:length(possible_subsets_of_paths)

		% %Run the Workers Asynchronously
		% for inner_subset_index = subset_index:subset_index+p.NumWorkers

		% end

		F{subset_index} = parfeval(p,@PerSubsetSynthesisAlgorithm,2,lcsas,subset_index,choices_as_binary_flags{subset_index},settings, possible_subsets_of_paths, MStar);

		if mod(subset_index,20) == 0
			disp(['Used parfeval ' num2str(subset_index) ' times.'])
		end

		% if optim0.problem == 0
		% 	break;
		% end

	end

	%Check to see if each activity has finished.
	numSubsetsFinished = 0;
	for subset_index = 1:length(F)
		numSubsetsFinished = numSubsetsFinished + strcmp(F{subset_index}.State,'finished');
	end

	next_update_at = 20;
	while numSubsetsFinished ~= length(possible_subsets_of_paths)

		numSubsetsFinished = 0;
		for subset_index = 1:length(F)
			if strcmp(F{subset_index}.State,'finished')
				numSubsetsFinished = numSubsetsFinished + 1;
			end
		end

		if mod(subset_index,20) == 0
			disp(['Used parfeval ' num2str(subset_index) ' times.'])
		end

		if numSubsetsFinished > next_update_at
			disp(['Completed ' num2str(numSubsetsFinished) ' attempts!' ])
			next_update_at = next_update_at + 20;
		else
			disp(['- No Updates. In Checking loop.'])
		end

	end

	% Fetch results
	for subset_index = 1:length(F)
		[ cbc , optim_info ] = fetchOutputs(F{subset_index});
		synthesis_info.Experiment{subset_index}.Controller = cbc;
		synthesis_info.Experiment{subset_index}.OptimizationInfo = optim_info;
	end


	synthesis_info.Timing.SearchForFeasibleP = toc(P_search_start);
	synthesis_info.Timing.Overall = toc(algorithm_start);

	controller = [];



end

function [ controller, optim0 ] = PerSubsetSynthesisAlgorithm( lcsas , subset_index , active_sequence_flags , settings , possible_subsets_of_paths , unchecked_knowledge_sequences )

	%% Constants

	L = lcsas.L;
	TimeHorizon = length(L.words{1});

	[ n_x , n_u , n_y , n_w , n_v ] = lcsas.Dimensions();
	cg = constr_gen(0);

	%% Algorithm

	%Start this iteration
	subset_index_start_time = tic;
	%active_sequence_flags = choices_as_binary_flags{subset_index};

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
								guaranteed_reachability_constraint + causal_detection_constraints + causal_gain_constraints;

	ops = sdpsettings(	'verbose',settings.verbosity,'debug',1);
	ops.gurobi.NodeLimit = 10^5;

	ops = sdpsettings(ops,'solver','gurobi');
	optim0 = optimize(optimization_constraints,[],ops);

	synthesis_info.Timing.Experiment{subset_index}.Optimization.yalmiptime = optim0.yalmiptime;
	synthesis_info.Timing.Experiment{subset_index}.Optimization.solvertime = optim0.solvertime; 
	synthesis_info.Timing.Experiment{subset_index}.ForExperiment = toc(subset_index_start_time);

end