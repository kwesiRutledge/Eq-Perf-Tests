function [ nonempty_flag , Kw_products , Hy_products , hy_products , constraints_out ] = get_nonempty_indicator_from_mccormick_envelopes(varargin)
	%Description:
	%	Creates the binary variable nonempty_flag which is one if the variable defined by H and h is nonempty
	%	and
	%
	%Usage:
	%	[ nonempty_flag , Hw_product , Hy_product , hy_product , constraints_out ] = cg.get_nonempty_indicator_from_mccormick_envelopes( lcsas_in , target_language_info , optimization_variables , grid_definition )
	%
	%Inputs:
	%	- target_language_info: A struct with the following fields: .L, .Internal_Behavior_H, .Internal_Behavior_h
	%	- optimization_variables: A struct containing the optimization variables for this problem (specifically for subproblem 3 as of right now)
	%							  Members include .K and .k
	%	- grid_definition: A struct containing the values that define the grids of interest for our McCormick Variables
	%					   Members include .NumberOfDivisionsPerKDimension, .NumberOfDivisionsPerWDimension, .NumberOfRegionsInH

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	[ cg , lcsas0 , target_language_info , optimization_variables , grid_definition , x0 , gnifme_settings ] = nonempty_indicator_from_mccormick_input_processing(varargin{:});

	
	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	L0 = target_language_info.L;

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	switch L0.cardinality() 
	 	case 1
	 		[ nonempty_flag , Kw_products , Hy_products , hy_products , constraints_out ] = get_nonempty_indicator_from_mccormick_envelopes_Card1(varargin{:});
	 	case 2
	 		[ nonempty_flag , Kw_products , Hy_products , hy_products , constraints_out ] = get_nonempty_indicator_from_mccormick_envelopes_Card2(varargin{:});
	 	otherwise
	 		error(['The function '])
	 end 

end

function [ nonempty_flag , Kw_products , Hy_products , hy_products , constraints_out ] = get_nonempty_indicator_from_mccormick_envelopes_Card1(varargin)
	%Description:
	%	Handles the version of this function when there is just one word in the target language.

	[ cg , lcsas0 , target_language_info , optimization_variables , grid_definition , x0 , gnifme_settings ] = nonempty_indicator_from_mccormick_input_processing(varargin{:});

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	K = optimization_variables.K;
	k = optimization_variables.k;

	H_IB = target_language_info.Internal_Behavior_H;
	h_IB = target_language_info.Internal_Behavior_h;

	L0 = target_language_info.L;

	TimeHorizon = length(L0.words{1});

	% The "word id" is the index of the words from L0 in terms of the indices of lcsas0.L
	%	For example, if lcsas0.L = { [1,2,3] , [4,5,6] , [7,8,9] } and
	%					L0       = { [1,2,3] , [7,8,9]}
	%	Then word_ids should be: [ 1 , 3 ] because the first word in L0 (i.e. L0.words{1}) is the first word of 
	%	lcsas0.L and the second word in L0 (i.e. L0.words{2}) is the THIRD word in lcsas0.L.

	% word_ids = [];
	% for word_index = 1:L0.cardinality()
	% 	[ ~ , temp_word_id ] = lcsas0.L.contains(L0.words{word_index});
	% 	word_ids(word_index) = temp_word_id;
	% end

	[ n_x , n_u , n_y , n_w , n_v ] = lcsas0.Dimensions();
	[ S_w , S_u , C_bar , J , f_bar ] = get_mpc_matrix_cell_arrays(lcsas0);

	% if length(word_ids) > 1
	% 	error('This function does not currently support languages with multiple words in it.')
	% end

	if L0.cardinality() > 1
		warning('This function has not been fully tested with languages with multiple words in it.')
	end

	NumberOfDivisionsPerKDimension = grid_definition.NumberOfDivisionsPerKDimension;
	NumberOfDivisionsPerkDimension = grid_definition.NumberOfDivisionsPerkDimension;
	NumberOfDivisionsPerWDimension = grid_definition.NumberOfDivisionsPerWDimension;
	NumberOfDivisionsPerHDimension = grid_definition.NumberOfDivisionsPerHDimension;
	NumberOfDivisionsPerYDimension = grid_definition.NumberOfDivisionsPerYDimension;

	Pw_arr = target_language_info.Disturbance_Set_Array;

    eps0 = 10^(-4);
    
	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	nonempty_flag = binvar(1,1,'full');

	% intersection_H = H_PhiI * ...
	% 	[ (S_w{word_id}+S_u{word_id}*K{knowl_seq_index}) ;
	% 		K{knowl_seq_index};
	% 		eye(n_w*TimeHorizon);
	% 		zeros(n_x,n_w*TimeHorizon) ];

	w = {}; y = {};
	Kw_products = {};
	Hy_products = {}; hy_products = {};

	interm_H = {}; interm_h = {};

	constraints_out = [];
	for word_index = 1:L0.cardinality()

		[~,word_id] = lcsas0.L.contains(L0.words{word_index});

		w{end+1} = sdpvar(n_w*TimeHorizon,1,'full');
		y{end+1} = sdpvar(size(H_IB,1),1,'full');

		%% Create the Matrices which Define "Allowable w's" Polytope

		nonK_prefactor = [ 	S_w{word_id} ;
							zeros(n_u*TimeHorizon,n_w*TimeHorizon) ;
							eye(n_w*TimeHorizon) ;
							zeros(n_x,n_w*TimeHorizon) ];

		K_prefactor = [ S_u{word_id} ;
						eye(n_u*TimeHorizon);
						zeros(n_w*TimeHorizon,n_u*TimeHorizon);
						zeros(n_x,n_u*TimeHorizon) ];

		H = H_IB * ( nonK_prefactor + K_prefactor * K );
		
		% h = h_IB - ...
		% 	H_IB * [ J{word_id}*x0 + S_w{word_id}*f_bar{word_id} + S_u{word_id}*k ;
		% 			k ;
		% 			zeros(n_w*TimeHorizon,1) ; 
		% 			x0 ];

		h_independent_factor = h_IB - ...
			H_IB * [ J{word_id}*x0 + S_w{word_id}*f_bar{word_id} ;
					zeros(n_u*TimeHorizon,1) ;
					zeros(n_w*TimeHorizon,1) ; 
					x0 ];

		h_dependent_factor = - ...
			H_IB * [ S_u{word_id} ;
					eye(n_u*TimeHorizon) ;
					zeros(n_w*TimeHorizon,n_u*TimeHorizon) ; 
					zeros(n_x,n_u*TimeHorizon) ];

		h = h_independent_factor + h_dependent_factor*k;

		%% Introduce Variables to map the constraints on H.

		interm_H{word_index} = sdpvar(size(H,1),size(H,2),'full');
		interm_h{word_index} = sdpvar(size(h,1),1,'full');

		bridge_constraints = [ interm_H{word_index} == H ] + [ interm_h{word_index} == h ];
		bound_constraints_on_interms = [ -gnifme_settings.H_bound <= interm_H{word_index} <= gnifme_settings.H_bound ] + [ -gnifme_settings.h_bound <= interm_h{word_index} <= gnifme_settings.h_bound ];

		constraints_out = constraints_out + bridge_constraints + bound_constraints_on_interms;

		%% Create McCormick Constraints for Products of Variables

		% Kw product constraints
		if gnifme_settings.debug_flag > 1
			disp('Kw product constraints')
		end

		grid_def_Kw = create_grid_def_for_Kw( lcsas0 , grid_definition , optimization_variables, target_language_info , gnifme_settings.K_bound , word_index );

		[ Kw_products{word_index} , mccormick_constraints , mccormick_binvars ] = cg.create_mccormick_envelope(	K,w{end},grid_def_Kw, ...
																												'debug_flag',gnifme_settings.debug_flag, ...
																												'eta_z_bounds',gnifme_settings.z_bound);
		constraints_out = constraints_out + mccormick_constraints;

		%Hy product constraints
		if gnifme_settings.debug_flag > 1
			disp('Hy product constraints')
		end

		grid_def_Hy = create_grid_def_for_Hy( lcsas0 , grid_definition , optimization_variables, target_language_info , gnifme_settings );

		[ Hy_products{word_index} , mccormick_constraints , mccormick_binvars ] = cg.create_mccormick_envelope(	interm_H{end}',y{end},grid_def_Hy, ...
																												'debug_flag',gnifme_settings.debug_flag, ...
																												'eta_z_bounds',gnifme_settings.z_bound );
		constraints_out = constraints_out + mccormick_constraints;

		%hy product constraints
		if gnifme_settings.debug_flag > 1
			disp('hy product constraints')
		end

		grid_def_hy = create_grid_def_for_hy( lcsas0 , grid_definition , optimization_variables , target_language_info , gnifme_settings );

		[ hy_products{word_index} , mccormick_constraints , mccormick_binvars ] = cg.create_mccormick_envelope(interm_h{end},y{end},grid_def_hy,...
																												'debug_flag',gnifme_settings.debug_flag, ...
																												'eta_z_bounds',gnifme_settings.z_bound );
		constraints_out = constraints_out + mccormick_constraints;

		%% Create the Containment Constraints

		constraints_out = constraints_out + ...
					[ Pw_arr(word_index).A*w{word_index} <= Pw_arr(word_index).b ] + ...
					implies( [nonempty_flag == 1] , [H_IB * ( nonK_prefactor*w{word_index} + K_prefactor * Kw_products{word_index} ) <= h]) + ...
					implies( [nonempty_flag == 0] , [Hy_products{word_index} == 0] + [ hy_products{word_index} <= -eps0 ] + [y{end} >= 0] );

	end

end

function [ nonempty_flag , Kw_products , Hy_products , hy_products , constraints_out ] = get_nonempty_indicator_from_mccormick_envelopes_Card2(varargin)
	%Description:
	%	Handles the version of this function when there is just one word in the target language.

	[ cg , lcsas0 , target_language_info , optimization_variables , grid_definition , x0 , gnifme_settings ] = nonempty_indicator_from_mccormick_input_processing(varargin{:});

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	K = optimization_variables.K;
	k = optimization_variables.k;

	H_IB = target_language_info.Internal_Behavior_H;
	h_IB = target_language_info.Internal_Behavior_h;

	L0 = target_language_info.L;
	card_L0 = L0.cardinality();

	TimeHorizon = length(L0.words{1});

	% The "word id" is the index of the words from L0 in terms of the indices of lcsas0.L
	%	For example, if lcsas0.L = { [1,2,3] , [4,5,6] , [7,8,9] } and
	%					L0       = { [1,2,3] , [7,8,9]}
	%	Then word_ids should be: [ 1 , 3 ] because the first word in L0 (i.e. L0.words{1}) is the first word of 
	%	lcsas0.L and the second word in L0 (i.e. L0.words{2}) is the THIRD word in lcsas0.L.

	% word_ids = [];
	% for word_index = 1:L0.cardinality()
	% 	[ ~ , temp_word_id ] = lcsas0.L.contains(L0.words{word_index});
	% 	word_ids(word_index) = temp_word_id;
	% end

	[ n_x , n_u , n_y , n_w , n_v ] = lcsas0.Dimensions();
	[ S_w , S_u , C_bar , J , f_bar ] = get_mpc_matrix_cell_arrays(lcsas0);

	% if length(word_ids) > 1
	% 	error('This function does not currently support languages with multiple words in it.')
	% end

	if L0.cardinality() > 1
		warning('This function has not been fully tested with languages with multiple words in it.')
	end

	NumberOfDivisionsPerKDimension = grid_definition.NumberOfDivisionsPerKDimension;
	NumberOfDivisionsPerkDimension = grid_definition.NumberOfDivisionsPerkDimension;
	NumberOfDivisionsPerWDimension = grid_definition.NumberOfDivisionsPerWDimension;
	NumberOfDivisionsPerHDimension = grid_definition.NumberOfDivisionsPerHDimension;
	NumberOfDivisionsPerYDimension = grid_definition.NumberOfDivisionsPerYDimension;

	Pw_arr = target_language_info.Disturbance_Set_Array;

    eps0 = 10^(-4);
    
	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	nonempty_flag = binvar(1,1,'full');

	% intersection_H = H_PhiI * ...
	% 	[ (S_w{word_id}+S_u{word_id}*K{knowl_seq_index}) ;
	% 		K{knowl_seq_index};
	% 		eye(n_w*TimeHorizon);
	% 		zeros(n_x,n_w*TimeHorizon) ];

	w = {}; y = {};
	Kw_products = {};
	Hy_products = {}; hy_products = {};

	interm_H = {}; interm_h = {};

	constraints_out = [];

	for word_index = 1:1

		[~,word_id] = lcsas0.L.contains(L0.words{1});

		w{end+1} = sdpvar(L0.cardinality()*n_w*TimeHorizon,1,'full');
		y{end+1} = sdpvar(size(H_IB,1),1,'full');

		%% Create the Matrices which Define "Allowable w's" Polytope

		nonK_prefactor1 = [ zeros(n_x*(TimeHorizon+1),(0)*n_w*TimeHorizon), S_w{word_id}, zeros(n_x*(TimeHorizon+1),(1)*n_w*TimeHorizon) ;
							zeros(n_u*TimeHorizon,card_L0*n_w*TimeHorizon) ;
							eye(card_L0*n_w*TimeHorizon) ;
							zeros(card_L0*n_x,card_L0*n_w*TimeHorizon) ];

		K_prefactor1 = [	S_u{word_id} ;
							eye(n_u*TimeHorizon);
							zeros(card_L0*n_w*TimeHorizon,n_u*TimeHorizon);
							zeros(card_L0*n_x,n_u*TimeHorizon) ] ;

		H1 = H_IB * ( nonK_prefactor1 + K_prefactor1 * K * [ eye(n_w*TimeHorizon) , zeros(n_w*TimeHorizon) ] ) ;
			

		% nonK_prefactor2 = [ zeros(n_w*TimeHorizon,(1)*n_w*TimeHorizon), S_w{word_id}, zeros((0)*n_w*TimeHorizon) ;
		% 					zeros(n_u*TimeHorizon,card_L0*n_w*TimeHorizon) ;
		% 					eye(card_L0*n_w*TimeHorizon) ;
		% 					zeros(card_L0*n_x,n_w*TimeHorizon) ];

		% K_prefactor2 = [	S_u{2} ;
		% 					eye(n_u*TimeHorizon);
		% 					zeros(card_L0*n_w*TimeHorizon,n_u*TimeHorizon);
		% 					zeros(card_L0*n_x,n_u*TimeHorizon) ] ;

		% H{2} = H_IB * ( nonK_prefactor + K_prefactor * K * [ zeros(n_w*TimeHorizon) , eye(n_w*TimeHorizon) ] ) ;

		% h = h_IB - ...
		% 	H_IB * [ J{word_id}*x0 + S_w{word_id}*f_bar{word_id} + S_u{word_id}*k ;
		% 			k ;
		% 			zeros(n_w*TimeHorizon,1) ; 
		% 			x0 ];

		h_independent_factor1 = h_IB - ...
			H_IB * [ J{word_id}*x0 + S_w{word_id}*f_bar{word_id} ;
					zeros(n_u*TimeHorizon,1) ;
					zeros(card_L0*n_w*TimeHorizon,1) ; 
					repmat(x0,card_L0,1) ];

		h_dependent_factor1 = - ...
			H_IB * [ S_u{word_id} ;
					eye(n_u*TimeHorizon) ;
					zeros(card_L0*n_w*TimeHorizon,n_u*TimeHorizon) ; 
					zeros(card_L0*n_x,n_u*TimeHorizon) ];

		h1 = h_independent_factor1 + h_dependent_factor1*k;

		%% Introduce Variables to map the constraints on H.

		interm_H{word_index} = sdpvar(size(H1,1),size(H1,2),'full');
		interm_h{word_index} = sdpvar(size(h1,1),1,'full');

		bridge_constraints = [ interm_H{word_index} == H1 ] + [ interm_h{word_index} == h1 ];
		constraints_out = constraints_out + bridge_constraints;

		%% Create McCormick Constraints for Products of Variables

		% Kw product constraints
		if gnifme_settings.debug_flag > 1
			disp('Kw product constraints')
		end

		grid_def_Kw = create_grid_def_for_Kw( lcsas0 , grid_definition , optimization_variables, target_language_info , gnifme_settings.K_bound , word_index );

		[ Kw_products{word_index} , mccormick_constraints , mccormick_binvars ] = cg.create_mccormick_envelope(	K,w{end}([1:n_w*TimeHorizon],1),grid_def_Kw,...
																												'debug_flag',gnifme_settings.debug_flag, ...
																												'eta_z_bounds',gnifme_settings.z_bound );
		constraints_out = constraints_out + mccormick_constraints;

		%Hy product constraints
		if gnifme_settings.debug_flag > 1
			disp('Hy product constraints')
		end

		grid_def_Hy = create_grid_def_for_Hy( lcsas0 , grid_definition , optimization_variables, target_language_info , gnifme_settings );

		[ Hy_products{word_index} , mccormick_constraints , mccormick_binvars ] = cg.create_mccormick_envelope(	interm_H{end}',y{end},grid_def_Hy, ...
																												'debug_flag',gnifme_settings.debug_flag, ...
																												'eta_z_bounds',gnifme_settings.z_bound );
		constraints_out = constraints_out + mccormick_constraints;

		%hy product constraints
		if gnifme_settings.debug_flag > 1
			disp('hy product constraints')
		end

		grid_def_hy = create_grid_def_for_hy( lcsas0 , grid_definition , optimization_variables , target_language_info , gnifme_settings );

		[ hy_products{word_index} , mccormick_constraints , mccormick_binvars ] = cg.create_mccormick_envelope( interm_h{end},y{end},grid_def_hy, ...
																												'debug_flag',gnifme_settings.debug_flag, ...
																												'eta_z_bounds',gnifme_settings.z_bound);
		constraints_out = constraints_out + mccormick_constraints;

		%% Create the Containment Constraints

		constraints_out = constraints_out + ...
					[ blkdiag(Pw_arr(1).A,Pw_arr(2).A)*w{word_index} <= [ Pw_arr(1).b ; Pw_arr(2).b ] ] + ...
					implies( [nonempty_flag == 1] , [H_IB * ( nonK_prefactor1*w{word_index} + K_prefactor1 * Kw_products{word_index} ) <= h1]) + ...
					implies( [nonempty_flag == 0] , [Hy_products{word_index} == 0] + [ hy_products{word_index} <= -eps0 ] + [y{end} >= 0] );

	end

end

function [ cg , lcsas_out , target_language_info , optimization_variables , grid_struct , x0 , gnifme_settings ] = nonempty_indicator_from_mccormick_input_processing(varargin)
	%Description:

	%% Constants %%

	eps1 = 10^5;
	target_language_info_fieldnames = {'L','Internal_Behavior_H','Internal_Behavior_h','Disturbance_Set_Array'};
	grid_struct_fieldnames = {'NumberOfDivisionsPerKDimension','NumberOfDivisionsPerWDimension','NumberOfDivisionsPerHDimension','NumberOfDivisionsPerYDimension'};

	%% Input Processing %%

	cg = varargin{1};
	lcsas_out = varargin{2};
	target_language_info = varargin{3};
	optimization_variables = varargin{4};
	grid_struct = varargin{5};
	x0 = varargin{6};

	%%%%%%%%%%%%%
	%% Default %%
	%%%%%%%%%%%%%

	gnifme_settings = struct('debug_flag',2,'K_bound',eps1,'H_bound',eps1,'y_bound',eps1,'h_bound',eps1,'z_bound',eps1);

	% grid_struct.H_lb = -eps1*ones(size(H));
	% grid_struct.H_ub =  eps1*ones(size(H));
	% grid_struct.h_lb = -eps1*ones(size(h));
	% grid_struct.h_ub =  eps1*ones(size(h));

	% grid_struct.w_lb = -eps1*ones(size(H,2),1);
	% grid_struct.w_ub =  eps1*ones(size(H,2),1);

	% grid_struct.y_lb = zero(size(H,1),1);
	% grid_struct.y_lb = eps1*ones(size(H,1),1);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Process Any Extra Inputs %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if nargin > 6
		argidx = 7;
		while argidx <= nargin
			switch varargin{argidx}
				case {'debug','debug_flag'}
					gnifme_settings.debug_flag = varargin{argidx+1};
					argidx = argidx + 2;
				case 'K_bound'
					gnifme_settings.K_bound = varargin{argidx+1};
					argidx = argidx + 2;
				case 'H_bound'
					gnifme_settings.H_bound = varargin{argidx+1};
					argidx = argidx + 2;
				case 'h_bound'
					gnifme_settings.h_bound = varargin{argidx+1};
					argidx = argidx + 2;
				case 'y_bound'
					gnifme_settings.y_bound = varargin{argidx+1};
					argidx = argidx + 2;
				otherwise
					error(['Unexpected input to get_nonempty_indicator_from_mccormick_envelopes: ' varargin{argidx} ])
			end
		end
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Test target_language_info %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%Verify that it has all of the desired fields.
	for fieldname_index = 1:length(target_language_info_fieldnames)

		if ~isfield( target_language_info , target_language_info_fieldnames{fieldname_index} )
			error([ 'The field name ' target_language_info_fieldnames{fieldname_index} ' is not part of the input grid_definition.' ])
		end

	end

	L0 = target_language_info.L;
	word_len1 = length(L0.words{1});

	for word_index = 2:L0.cardinality()
		temp_word = L0.words{word_index};
		if length(temp_word) ~= word_len1
			error(['This function supports languages that are all of the same length. Word 1 is of length ' num2str(word_len1) ' while word ' num2str(word_index) ' is of length ' num2str(length(temp_word)) '.' ])
		end
	end

	%%%%%%%%%%%%%%%%%%%%%%
	%% Test grid_struct %%
	%%%%%%%%%%%%%%%%%%%%%%

	%Verify that it has all of the desired fields.
	for fieldname_index = 1:length(grid_struct_fieldnames)

		if ~isfield( grid_struct , grid_struct_fieldnames{fieldname_index} )
			error([ 'The field name ' grid_struct_fieldnames{fieldname_index} ' is not part of the input grid_definition.' ])
		end

	end


end

function [ S_w , S_u , C_bar , J , f_bar ] = get_mpc_matrix_cell_arrays(lcsas_in)

	%% Constants

	card_L = lcsas_in.L.cardinality();

	%% Algorithm
	S_w = { };
	S_u = { };
	C_bar = { };
	J = { };
	f_bar = { };

	for word_index = 1:card_L
		[S_w{word_index},S_u{word_index},C_bar{word_index},J{word_index},f_bar{word_index}] = get_mpc_matrices(lcsas_in,'word',lcsas_in.L.words{word_index});
	end

end

function [ grid_def_out ] = create_grid_def_for_Kw( lcsas_in , grid_struct , optimization_variables , target_language_info , K_bound , current_word_index )
	%Description:
	%	Creates a grid definition for the Kw product which takes into account that:
	%	- K is lower block triangular
	%	- 

	%% Constants

	[ n_x , n_u , n_y , n_w , n_v ] = lcsas_in.Dimensions();
	L0 = target_language_info.L;
	T = length(L0.words{1});

	Pw_arr = target_language_info.Disturbance_Set_Array;

	NumberOfDivisionsPerKDimension = grid_struct.NumberOfDivisionsPerKDimension;
	NumberOfDivisionsPerWDimension = grid_struct.NumberOfDivisionsPerWDimension;

	K = optimization_variables.K;

	%% Algorithm

	K_lb = -K_bound*ones(size(K));
	K_ub = K_bound*ones( size(K) );
	num_K_regions = NumberOfDivisionsPerKDimension*ones(size(K));
	for t = 1:T-1
		K_lb((t-1)*n_u+[1:n_u],[t*n_w+1:end]) = 0;
		K_ub((t-1)*n_u+[1:n_u],[t*n_w+1:end]) = 0;
		num_K_regions((t-1)*n_u+[1:n_u],[t*n_w+1:end]) = 1;
	end

	% Create lb for w

	grid_def_out = struct( ...
					'x_lb' , K_lb , 'x_ub' , K_ub , ...
					'y_lb' , Pw_arr(current_word_index).Internal.lb , 'y_ub' , Pw_arr(current_word_index).Internal.ub , ...
					'NumberOfRegionsInX', num_K_regions, ...
					'NumberOfRegionsInY', NumberOfDivisionsPerWDimension * ones(T*n_w));

end

function [ grid_def_out ] = create_grid_def_for_Hy( lcsas_in , grid_struct , optimization_variables , target_language_info , settings_struct )
	%Description:
	%	Creates a grid definition for the Kw product which takes into account that:
	%	- K is lower block triangular
	%	- 

	%% Constants

	[ n_x , n_u , n_y , n_w , n_v ] = lcsas_in.Dimensions();
	L0 = target_language_info.L;
	card_L0 = L0.cardinality();
	T = length(L0.words{1});

	NumberOfDivisionsPerHDimension = grid_struct.NumberOfDivisionsPerHDimension;
	NumberOfDivisionsPerYDimension = grid_struct.NumberOfDivisionsPerYDimension;

	K = optimization_variables.K;

	H_bound = settings_struct.H_bound;
	y_bound = settings_struct.y_bound;

	%% Algorithm

	H_lb = -H_bound*ones( size( target_language_info.Internal_Behavior_H,1 ) , card_L0*n_w*T )';
	H_ub = H_bound*ones( size( target_language_info.Internal_Behavior_H,1 ) , card_L0*n_w*T )';

	% Create lb for w

	grid_def_out = struct( ...
						'x_lb',H_lb,'x_ub',H_ub, ...
						'y_lb',zeros(size(H_lb,2),1),'y_ub',y_bound*ones(size(H_lb,2),1) , ...
						'NumberOfRegionsInX', NumberOfDivisionsPerHDimension*ones(size(H_lb)), ...
						'NumberOfRegionsInY', NumberOfDivisionsPerYDimension * ones(size(H_lb,2)));

end

function [ grid_def_out ] = create_grid_def_for_hy( lcsas_in , grid_struct , optimization_variables , target_language_info , settings_struct )
	%Description:
	%	Creates a grid definition for the Kw product which takes into account that:
	%	- K is lower block triangular
	%	- 

	%% Constants

	[ n_x , n_u , n_y , n_w , n_v ] = lcsas_in.Dimensions();
	L0 = target_language_info.L;
	card_L0 = L0.cardinality();
	T = length(L0.words{1});

	NumberOfDivisionsPerhDimension = grid_struct.NumberOfDivisionsPerhDimension;
	NumberOfDivisionsPerYDimension = grid_struct.NumberOfDivisionsPerYDimension;

	h_bound = settings_struct.h_bound;
	y_bound = settings_struct.y_bound;

	%% Algorithm

	h_lb = -h_bound*ones( size( target_language_info.Internal_Behavior_h,1 ) , 1 );
	h_ub = h_bound*ones( size( target_language_info.Internal_Behavior_h,1 ) , 1 );

	y_ub = y_bound*ones(size(h_lb,1),1);

	% Create lb for w

	grid_def_out = struct( ...
						'x_lb',h_lb,'x_ub',h_ub, ...
						'y_lb',zeros(size(h_lb,1),1),'y_ub',y_ub , ...
						'NumberOfRegionsInX', NumberOfDivisionsPerhDimension*ones(size(h_lb,1)), ...
						'NumberOfRegionsInY', NumberOfDivisionsPerYDimension * ones(size(h_lb,1)));

end