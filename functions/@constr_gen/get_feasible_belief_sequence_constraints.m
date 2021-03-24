function [ varargout ] = get_feasible_belief_sequence_constraints( varargin )
	%Description:
	%	Creates constraints which represent whether or not the the knowledge sequences in knowl_seq_matrix
	%	are feasible given the choice of K_cell_arr , k_cell_arr and the selected assignment of empty (0)
	%	or non-empty (1) given in selected_knowl_sequences
	%
	%Usage:
	%	[ feasibility_constraints , feas_constr_vars ] = cg.get_feasible_belief_sequence_constraints( lcsas0 , knowl_seq_matrix , K_cell_arr , k_cell_arr , selected_knowl_sequences , Pu , ConsistentDisturbs , x0 )
	%	[ feasibility_constraints , feas_constr_vars ] = cg.get_feasible_belief_sequence_constraints( lcsas0 , knowl_seq_matrix , K_cell_arr , k_cell_arr , selected_knowl_sequences , Pu , ConsistentDisturbs , x0 , 'ConstraintFormat' , 'Gurobi')

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	[ ~ , ~ , ~ , ~ , ~ , ~ , ~ , ~ , ~ , settings0 ] = input_processing_gfbsc(varargin{:});

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	switch settings0.ConstraintFormat
		case 'YALMIP'
			[ varargout{1} , varargout{2} ] = get_feasible_belief_sequence_constraints_in_YALMIP( varargin{:} );
		case 'Gurobi'
			[ varargout{1} , varargout{2} ] = get_feasible_belief_sequence_constraints_in_Gurobi( varargin{:} );
			
		otherwise
			error(['Unexpected constraint format given: ' settings0.ConstraintFormat ])
	end
		

end

function [ cg , lcsas0 , knowl_seq_matrix , K_cell_arr , k_cell_arr , selected_knowl_sequences , Pu , ConsistentDisturbs, x0 , settings0 ] = input_processing_gfbsc(varargin)
	%Description:
	%	Parsing through the information available in the inputs to get_feasible_belief_sequence_constraints() and doing some sanity checks.


	% Check the number of inputs
	if nargin < 9
		error(['Not enough input arguments. Expecting at least 7; received ' num2str(nargin) '.' ])
	end

	%%%%%%%%%%%%%%%%%%%%%
	%% Retrieve values %%
	%%%%%%%%%%%%%%%%%%%%%

	cg = varargin{1};
	lcsas0 = varargin{2};
	knowl_seq_matrix = varargin{3};
	K_cell_arr = varargin{4};
	k_cell_arr = varargin{5};
	selected_knowl_sequences = varargin{6};
	Pu = varargin{7};
	ConsistentDisturbs = varargin{8};
	x0 = varargin{9};

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Check for any additional components %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	settings0 = [];

	if nargin > 9
		argidx = 10;
		while argidx <= nargin
			switch varargin{argidx}
				case 'ConstraintFormat'
					settings0.ConstraintFormat = varargin{argidx+1};
					argidx = argidx + 2;
				otherwise
					error(['Unexpected input to function: ' varargin{argidx} ])
			end
		end
	end


	%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Default Definitions %%
	%%%%%%%%%%%%%%%%%%%%%%%%%

	if ~isfield(settings0,'ConstraintFormat')
		settings0.ConstraintFormat = 'YALMIP';
	end

	%%%%%%%%%%%
	%% LCSAS %%
	%%%%%%%%%%%

	% Check that LCSAS contains nonempty X0 value.
	if isempty(lcsas0.X0)
		error('The value of lcsas0.X0 must be nonempty. Please assign it to be a Polyhedron element before calling this function.')
	end

	%%%%%%%%
	%% Pu %%
	%%%%%%%%

	if ~isa(Pu,'Polyhedron')
		error(['Expected Pu to be of class Polyhedron; instead it is of class ' class(Pu) '.'])
	end

end

function [ constraints_out , new_variables_out ] = get_feasible_belief_sequence_constraints_in_YALMIP( varargin )
	%Description:
	%	Creates constraints which represent whether or not the the knowledge sequences in knowl_seq_matrix
	%	are feasible given the choice of K_cell_arr , k_cell_arr and the selected assignment of empty (0)
	%	or non-empty (1) given in selected_knowl_sequences
	%
	%Usage:
	%	[ feasibility_constraints , feas_constr_vars ] = cg.get_feasible_belief_sequence_constraints( lcsas0 , knowl_seq_matrix , K_cell_arr , k_cell_arr , selected_knowl_sequences , Pu , ConsistentDisturbs , x0 )

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	[ cg , lcsas0 , knowl_seq_matrix , K_cell_arr , k_cell_arr , selected_knowl_sequences , Pu , ConsistentDisturbs , x0 ] = input_processing_gfbsc(varargin{:});

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	TimeHorizon = size(knowl_seq_matrix,1);
	num_knowl_sequences = size(knowl_seq_matrix,2);

	[ n_x , n_u , n_y , n_w , n_v ] = lcsas0.Dimensions();

	[ S_w , S_u , ~ , J , f_bar ] = lcsas0.get_mpc_matrices('All Words');

	eps0 = 10^(-5);

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	matching_behavior_constraint = [];
	dummy_var_bound_constraints = []; temp_dummy_w = {}; temp_dummy_y = {};
	for knowl_seq_index = 1:num_knowl_sequences
		
		knowl_seq = knowl_seq_matrix(:,knowl_seq_index);
		temp_last_lang = knowl_seq(end);
		tll_card = temp_last_lang.cardinality();

		%Get internal_behavior_matrices
		[ H_PhiI , h_PhiI ] = lcsas0.get_consistent_internal_behavior_matrices( TimeHorizon , knowl_seq(end) , Pu , lcsas0.X0 );

		nonK_prefactor = {}; K_prefactor = {}; h_independent_factor = {}; h_dependent_factor = {};
		H = {}; h = {};
		switch temp_last_lang.cardinality()
		 	case 1
		 	% 	intersection_H = [ H_PhiI ];
				% intersection_H = [ 	intersection_H ;
				% 					zeros( n_u*t ,n_x*(t+1)) , -eye(n_u*t) , K{gain_index}([1:n_u*t],[1:n_w*t]) , zeros(n_u*t,n_x) ;
				% 					zeros( n_u*t ,n_x*(t+1)) , eye(n_u*t) , -K{gain_index}([1:n_u*t],[1:n_w*t]) , zeros(n_u*t,n_x) ];

				[~,word_id] = lcsas0.L.contains(temp_last_lang.words{1});

				nonK_prefactor{1} = [ 	S_w{word_id} ;
					zeros(n_u*TimeHorizon,n_w*TimeHorizon) ;
					eye(n_w*TimeHorizon) ;
					zeros(n_x,n_w*TimeHorizon) ];

				K_prefactor{1} = [ S_u{word_id} ;
								eye(n_u*TimeHorizon);
								zeros(n_w*TimeHorizon,n_u*TimeHorizon);
								zeros(n_x,n_u*TimeHorizon) ];

				H{1} = H_PhiI * ( nonK_prefactor{1} + K_prefactor{1}*K_cell_arr{knowl_seq_index} );

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

				h{1} = h_independent_factor{1} + h_dependent_factor{1}*k_cell_arr{knowl_seq_index};

				temp_dummy_w{end+1} = sdpvar(size(H{1},2),1,'full');
				temp_dummy_y{end+1} = sdpvar(size(H{1},1),1,'full');

				% matching_behavior_constraint = matching_behavior_constraint + [0 <= matching_behavior{t+1}(knowl_seq_index) <= 1] + ...
				% 	iff(intersection_H * temp_dummy_eta{end} <= intersection_h, matching_behavior{t+1}(knowl_seq_index) == 1) + ...
				% 	iff([intersection_H'*temp_dummy_y{end} == 0] + [ intersection_h'*temp_dummy_y{end} <= -eps0 ] + [temp_dummy_y{end} >= 0] , matching_behavior{t+1}(knowl_seq_index) == 0);
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

					H{tll_index} = H_PhiI * ( nonK_prefactor{tll_index} + ...
						K_prefactor{tll_index} * K_cell_arr{knowl_seq_index} * [ zeros(n_w*TimeHorizon,(tll_index-1)*n_w*TimeHorizon), eye(n_w*TimeHorizon), zeros(n_w*TimeHorizon,(tll_card-tll_index)*n_w*TimeHorizon) ] ) ;
					
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

					h{tll_index} = h_independent_factor{tll_index} + h_dependent_factor{tll_index}*k_cell_arr{knowl_seq_index};

					temp_dummy_w{end+1} = sdpvar(size(H{tll_index},2),1,'full');
					temp_dummy_y{end+1} = sdpvar(size(H{tll_index},1),1,'full');

					% % Create Constraints
					% dummy_var_bound_constraints = dummy_var_bound_constraints + [ -temp_K_bound <= temp_dummy_eta{end} <= temp_K_bound ] + [ temp_dummy_y{end} <= temp_K_bound ];

				end

				% matching_behavior_constraint = matching_behavior_constraint + [0 <= matching_behavior{t+1}(knowl_seq_index) <= 1] + ...
				% 	iff(intersection_H * temp_dummy_eta{end} <= intersection_h, matching_behavior{t+1}(knowl_seq_index) == 1) + ...
				% 	iff([intersection_H'*temp_dummy_y{end} == 0] + [ intersection_h'*temp_dummy_y{end} <= -eps0 ] + [temp_dummy_y{end} >= 0] , matching_behavior{t+1}(knowl_seq_index) == 0);
		 	otherwise
		 		error("Unexpected cardinality.")
		 end 

		 % Create Constraints
		 if selected_knowl_sequences(knowl_seq_index)
		 	Pw_prime = 1;
		 	for Pw_index = 1:size(ConsistentDisturbs{knowl_seq_index})
		 		Pw_prime = Pw_prime * ConsistentDisturbs{knowl_seq_index}(Pw_index);
		 	end

			%Apply the existence of w condition.
			matching_behavior_constraint = matching_behavior_constraint + [ H{1} * temp_dummy_w{end} <= h{1} ] + [ Pw_prime.A * temp_dummy_w{end} <= Pw_prime.b ];
		else
			%Apply the empty polyhedron constraint
			matching_behavior_constraint = matching_behavior_constraint + ...
				[H{1}'*temp_dummy_y{end} == 0] + [ h{1}'*temp_dummy_y{end} <= -eps0 ] + [temp_dummy_y{end} >= 0];
		end

	end

	constraints_out = matching_behavior_constraint;
	new_variables_out = {temp_dummy_w,temp_dummy_y};

	% Create Outputs
	

end

function [ constraints_out , new_variables_out ] = get_feasible_belief_sequence_constraints_in_Gurobi( varargin )
	%Description:
	%	Creates constraints which represent whether or not the the knowledge sequences in knowl_seq_matrix
	%	are feasible given the choice of K_cell_arr , k_cell_arr and the selected assignment of empty (0)
	%	or non-empty (1) given in selected_knowl_sequences
	%
	%Usage:
	%	[ feasibility_constraints , feas_constr_vars ] = cg.get_feasible_belief_sequence_constraints( lcsas0 , knowl_seq_matrix , K_cell_arr , k_cell_arr , selected_knowl_sequences , Pu , ConsistentDisturbs , x0 )

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	[ cg , lcsas0 , knowl_seq_matrix , K_cell_arr , k_cell_arr , selected_knowl_sequences , Pu , ConsistentDisturbs , x0 ] = input_processing_gfbsc(varargin{:});

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	TimeHorizon = size(knowl_seq_matrix,1);
	num_knowl_sequences = size(knowl_seq_matrix,2);

	[ n_x , n_u , n_y , n_w , n_v ] = lcsas0.Dimensions();

	[ S_w , S_u , ~ , J , f_bar ] = lcsas0.get_mpc_matrices('All Words');

	eps0 = 10^(-5);

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	matching_behavior_constraint = [];
	dummy_var_bound_constraints = []; temp_dummy_w = {}; temp_dummy_y = {};
	for knowl_seq_index = 1:num_knowl_sequences
		
		knowl_seq = knowl_seq_matrix(:,knowl_seq_index);
		temp_last_lang = knowl_seq(end);
		tll_card = temp_last_lang.cardinality();

		%Get internal_behavior_matrices
		[ H_PhiI , h_PhiI ] = lcsas0.get_consistent_internal_behavior_matrices( TimeHorizon , knowl_seq(end) , Pu , lcsas0.X0 );

		nonK_prefactor = {}; K_prefactor = {}; h_independent_factor = {}; h_dependent_factor = {};
		H = {}; h = {};
		switch temp_last_lang.cardinality()
		 	case 1
		 	% 	intersection_H = [ H_PhiI ];
				% intersection_H = [ 	intersection_H ;
				% 					zeros( n_u*t ,n_x*(t+1)) , -eye(n_u*t) , K{gain_index}([1:n_u*t],[1:n_w*t]) , zeros(n_u*t,n_x) ;
				% 					zeros( n_u*t ,n_x*(t+1)) , eye(n_u*t) , -K{gain_index}([1:n_u*t],[1:n_w*t]) , zeros(n_u*t,n_x) ];

				[~,word_id] = lcsas0.L.contains(temp_last_lang.words{1});

				nonK_prefactor{1} = [ 	S_w{word_id} ;
					zeros(n_u*TimeHorizon,n_w*TimeHorizon) ;
					eye(n_w*TimeHorizon) ;
					zeros(n_x,n_w*TimeHorizon) ];

				K_prefactor{1} = [ S_u{word_id} ;
								eye(n_u*TimeHorizon);
								zeros(n_w*TimeHorizon,n_u*TimeHorizon);
								zeros(n_x,n_u*TimeHorizon) ];

				H{1} = H_PhiI * ( nonK_prefactor{1} + K_prefactor{1}*K_cell_arr{knowl_seq_index} );

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

				h{1} = h_independent_factor{1} + h_dependent_factor{1}*k_cell_arr{knowl_seq_index};

				temp_dummy_w{end+1} = sdpvar(size(H{1},2),1,'full');
				temp_dummy_y{end+1} = sdpvar(size(H{1},1),1,'full');

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

					H{tll_index} = H_PhiI * ( nonK_prefactor{tll_index} + ...
						K_prefactor{tll_index} * K_cell_arr{knowl_seq_index} * [ zeros(n_w*TimeHorizon,(tll_index-1)*n_w*TimeHorizon), eye(n_w*TimeHorizon), zeros(n_w*TimeHorizon,(tll_card-tll_index)*n_w*TimeHorizon) ] ) ;
					
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

					h{tll_index} = h_independent_factor{tll_index} + h_dependent_factor{tll_index}*k_cell_arr{knowl_seq_index};

					temp_dummy_w{end+1} = sdpvar(size(H{tll_index},2),1,'full');
					temp_dummy_y{end+1} = sdpvar(size(H{tll_index},1),1,'full');

					% % Create Constraints
					% dummy_var_bound_constraints = dummy_var_bound_constraints + [ -temp_K_bound <= temp_dummy_eta{end} <= temp_K_bound ] + [ temp_dummy_y{end} <= temp_K_bound ];

				end

				% matching_behavior_constraint = matching_behavior_constraint + [0 <= matching_behavior{t+1}(knowl_seq_index) <= 1] + ...
				% 	iff(intersection_H * temp_dummy_eta{end} <= intersection_h, matching_behavior{t+1}(knowl_seq_index) == 1) + ...
				% 	iff([intersection_H'*temp_dummy_y{end} == 0] + [ intersection_h'*temp_dummy_y{end} <= -eps0 ] + [temp_dummy_y{end} >= 0] , matching_behavior{t+1}(knowl_seq_index) == 0);
		 	otherwise
		 		error("Unexpected cardinality.")
		 end 

		 % Create Constraints
		 if selected_knowl_sequences(knowl_seq_index)
		 	Pw_prime = 1;
		 	for Pw_index = 1:size(ConsistentDisturbs{knowl_seq_index})
		 		Pw_prime = Pw_prime * ConsistentDisturbs{knowl_seq_index}(Pw_index);
		 	end

			%Apply the existence of w condition.
			matching_behavior_constraint = matching_behavior_constraint + [ H{1} * temp_dummy_w{end} <= h{1} ] + [ Pw_prime.A * temp_dummy_w{end} <= Pw_prime.b ];
		else
			%Apply the empty polyhedron constraint
			matching_behavior_constraint = matching_behavior_constraint + ...
				[H{1}'*temp_dummy_y{end} == 0] + [ h{1}'*temp_dummy_y{end} <= -eps0 ] + [temp_dummy_y{end} >= 0];
		end

	end

	constraints_out = matching_behavior_constraint;
	new_variables_out = {temp_dummy_w,temp_dummy_y};

end

function [ constraints_out , new_variables_out ] = get_feasible_belief_sequence_constraints_in_Symbolic( varargin )
	%Description:
	%	Creates constraints which represent whether or not the the knowledge sequences in knowl_seq_matrix
	%	are feasible given the choice of K_cell_arr , k_cell_arr and the selected assignment of empty (0)
	%	or non-empty (1) given in selected_knowl_sequences
	%
	%Usage:
	%	[ feasibility_constraints , feas_constr_vars ] = cg.get_feasible_belief_sequence_constraints( lcsas0 , knowl_seq_matrix , K_cell_arr , k_cell_arr , selected_knowl_sequences , Pu , ConsistentDisturbs , x0 )

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	[ cg , lcsas0 , knowl_seq_matrix , K_cell_arr , k_cell_arr , selected_knowl_sequences , Pu , ConsistentDisturbs , x0 ] = input_processing_gfbsc(varargin{:});

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	TimeHorizon = size(knowl_seq_matrix,1);
	num_knowl_sequences = size(knowl_seq_matrix,2);

	[ n_x , n_u , n_y , n_w , n_v ] = lcsas0.Dimensions();

	[ S_w , S_u , ~ , J , f_bar ] = lcsas0.get_mpc_matrices('All Words');

	eps0 = 10^(-5);

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	matching_behavior_constraint = [];
	dummy_var_bound_constraints = []; temp_dummy_w = {}; temp_dummy_y = {};
	for knowl_seq_index = 1:num_knowl_sequences
		
		knowl_seq = knowl_seq_matrix(:,knowl_seq_index);
		temp_last_lang = knowl_seq(end);
		tll_card = temp_last_lang.cardinality();

		%Get internal_behavior_matrices
		[ H_PhiI , h_PhiI ] = lcsas0.get_consistent_internal_behavior_matrices( TimeHorizon , knowl_seq(end) , Pu , lcsas0.X0 );

		nonK_prefactor = {}; K_prefactor = {}; h_independent_factor = {}; h_dependent_factor = {};
		H = {}; h = {};
		switch temp_last_lang.cardinality()
		 	case 1
		 	% 	intersection_H = [ H_PhiI ];
				% intersection_H = [ 	intersection_H ;
				% 					zeros( n_u*t ,n_x*(t+1)) , -eye(n_u*t) , K{gain_index}([1:n_u*t],[1:n_w*t]) , zeros(n_u*t,n_x) ;
				% 					zeros( n_u*t ,n_x*(t+1)) , eye(n_u*t) , -K{gain_index}([1:n_u*t],[1:n_w*t]) , zeros(n_u*t,n_x) ];

				[~,word_id] = lcsas0.L.contains(temp_last_lang.words{1});

				nonK_prefactor{1} = [ 	S_w{word_id} ;
					zeros(n_u*TimeHorizon,n_w*TimeHorizon) ;
					eye(n_w*TimeHorizon) ;
					zeros(n_x,n_w*TimeHorizon) ];

				K_prefactor{1} = [ S_u{word_id} ;
								eye(n_u*TimeHorizon);
								zeros(n_w*TimeHorizon,n_u*TimeHorizon);
								zeros(n_x,n_u*TimeHorizon) ];

				H{1} = H_PhiI * ( nonK_prefactor{1} + K_prefactor{1}*K_cell_arr{knowl_seq_index} );

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

				h{1} = h_independent_factor{1} + h_dependent_factor{1}*k_cell_arr{knowl_seq_index};

				temp_dummy_w{end+1} = sdpvar(size(H{1},2),1,'full');
				temp_dummy_y{end+1} = sdpvar(size(H{1},1),1,'full');

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

					H{tll_index} = H_PhiI * ( nonK_prefactor{tll_index} + ...
						K_prefactor{tll_index} * K_cell_arr{knowl_seq_index} * [ zeros(n_w*TimeHorizon,(tll_index-1)*n_w*TimeHorizon), eye(n_w*TimeHorizon), zeros(n_w*TimeHorizon,(tll_card-tll_index)*n_w*TimeHorizon) ] ) ;
					
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

					h{tll_index} = h_independent_factor{tll_index} + h_dependent_factor{tll_index}*k_cell_arr{knowl_seq_index};

					temp_dummy_w{end+1} = sdpvar(size(H{tll_index},2),1,'full');
					temp_dummy_y{end+1} = sdpvar(size(H{tll_index},1),1,'full');

					% % Create Constraints
					% dummy_var_bound_constraints = dummy_var_bound_constraints + [ -temp_K_bound <= temp_dummy_eta{end} <= temp_K_bound ] + [ temp_dummy_y{end} <= temp_K_bound ];

				end

				% matching_behavior_constraint = matching_behavior_constraint + [0 <= matching_behavior{t+1}(knowl_seq_index) <= 1] + ...
				% 	iff(intersection_H * temp_dummy_eta{end} <= intersection_h, matching_behavior{t+1}(knowl_seq_index) == 1) + ...
				% 	iff([intersection_H'*temp_dummy_y{end} == 0] + [ intersection_h'*temp_dummy_y{end} <= -eps0 ] + [temp_dummy_y{end} >= 0] , matching_behavior{t+1}(knowl_seq_index) == 0);
		 	otherwise
		 		error("Unexpected cardinality.")
		 end 

		 % Create Constraints
		 if selected_knowl_sequences(knowl_seq_index)
		 	Pw_prime = 1;
		 	for Pw_index = 1:size(ConsistentDisturbs{knowl_seq_index})
		 		Pw_prime = Pw_prime * ConsistentDisturbs{knowl_seq_index}(Pw_index);
		 	end

			%Apply the existence of w condition.
			matching_behavior_constraint = matching_behavior_constraint + [ H{1} * temp_dummy_w{end} <= h{1} ] + [ Pw_prime.A * temp_dummy_w{end} <= Pw_prime.b ];
		else
			%Apply the empty polyhedron constraint
			matching_behavior_constraint = matching_behavior_constraint + ...
				[H{1}'*temp_dummy_y{end} == 0] + [ h{1}'*temp_dummy_y{end} <= -eps0 ] + [temp_dummy_y{end} >= 0];
		end

	end

	constraints_out = matching_behavior_constraint;
	new_variables_out = {temp_dummy_w,temp_dummy_y};

end