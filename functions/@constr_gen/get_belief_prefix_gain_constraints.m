function [ constraints_out ] = get_belief_prefix_gain_constraints( varargin )
	%Description:
	%	
	%Usage:
	%	constrs = cg.get_belief_prefix_gain_constraints( lcsas0 , K_cell_arr , k_cell_arr , LK_sequences_in )
	%	constrs = cg.get_belief_prefix_gain_constraints( lcsas0 , K_cell_arr , k_cell_arr , LK_sequences_in , 'Ignore K')
	%	constrs = cg.get_belief_prefix_gain_constraints( lcsas0 , K_cell_arr , k_cell_arr , LK_sequences_in , 'Ignore k')

	%% Input Processing

	[ constraint_generator , lcsas_in , K_cell_arr , k_cell_arr , LK_sequences_in , gbpgc_settings ] = ip_get_belief_prefix_gain_constraints(varargin{:});

	%% Constants
	num_sequences = size(LK_sequences_in,2);

	[ n_x , n_u , n_y , n_w , n_v ] = lcsas_in.Dimensions();

	%% Algorithm

	% Determine Which Prefixes Match One Another
	shared_prefix_length = NaN(num_sequences) ;

	for knowl_seq_index = 1:num_sequences
		temp_knowl_seq = LK_sequences_in(:,knowl_seq_index);
		for knowl_seq_index2 = knowl_seq_index+1:num_sequences
			temp_knowl_seq2 = LK_sequences_in(:,knowl_seq_index2);
			
			shared_prefix_length(knowl_seq_index,knowl_seq_index2) = 1; %All pairs share at least one with one another.
			for prefix_length = 1:length(temp_knowl_seq2)
				if constraint_generator.knowledge_sequence_contains_prefix( temp_knowl_seq , temp_knowl_seq2([1:prefix_length],1) )
					shared_prefix_length(knowl_seq_index,knowl_seq_index2) = prefix_length;
				end
			end
		end
	end

	% Apply Constraints to the big gain matrix K_cell_arr
	constraints_out = [];

	if ~gbpgc_settings.ignore_K

		for knowl_seq_index = 1:num_sequences
			for knowl_seq_index2 = knowl_seq_index+1:num_sequences
				%Iterate through every row of the gain at K{knowl_seq_index}
				for prefix_length = 1:shared_prefix_length(knowl_seq_index,knowl_seq_index2)
					constraints_out = constraints_out + ...
						[ K_cell_arr{knowl_seq_index}((prefix_length-1)*n_u+[1:n_u],:) == K_cell_arr{knowl_seq_index2}((prefix_length-1)*n_u+[1:n_u],:) ];
				end
			end
		end

	end

	% Apply Constraints to the small gain matrix k_cell_arr

	if ~gbpgc_settings.ignore_k

		for knowl_seq_index = 1:num_sequences
			for knowl_seq_index2 = knowl_seq_index+1:num_sequences
				%Iterate through every row of the gain at K{knowl_seq_index}
				shared_prefix_length12 = shared_prefix_length(knowl_seq_index,knowl_seq_index2);
				constraints_out = constraints_out + ...
						[ k_cell_arr{knowl_seq_index}([1:(prefix_length)*n_u]) == k_cell_arr{knowl_seq_index2}([1:(prefix_length)*n_u]) ];
				
			end
		end
		
	end

end

function [ constraint_generator , lcsas_in , K_cell_arr , k_cell_arr , LK_sequences_in , gbpgc_settings ] = ip_get_belief_prefix_gain_constraints(varargin)
	%Description:
	%	Processes the many inputs to get_belief_prefix_gain_constraints().
	%
	%Usage:
	%	[ constraint_generator , lcsas_in , K_cell_arr , k_cell_arr , LK_sequences_in , gbpgc_settings ] = ip_get_belief_prefix_gain_constraints(varargin{:});

	%% Check If Enough Arguments are Provided
	if nargin < 5
		error(['get_belief_prefix_gain_constraints() expects at least 5 arguments, received ' num2str(nargin) ])
	end

	%% Collect First 5
	constraint_generator = varargin{1};
	lcsas_in = varargin{2};
	K_cell_arr = varargin{3};
	k_cell_arr = varargin{4};
	LK_sequences_in = varargin{5};

	%% Checking The Important 5
	if length(K_cell_arr) ~= size(LK_sequences_in,2)
		error('There should be the same number of columns in LK_sequences_in as in K_cell_arr.')
	end

	%% Collect Extras
	gbpgc_settings = struct( ...
		'ignore_k',false, ...
		'ignore_K',false ...
		);

	argument_index = 6;
	while argument_index <= nargin
		switch varargin{argument_index}
		case 'Ignore K'
			gbpgc_settings.ignore_K = true;
			argument_index = argument_index + 1;
		case 'Ignore k'
			gbpgc_settings.ignore_k = true;
			argument_index = argument_index + 1;
		otherwise
			error(['Unexpected input to get_belief_prefix_gain_constraints(): ' varargin{argument_index} ])
		end 
	end

end