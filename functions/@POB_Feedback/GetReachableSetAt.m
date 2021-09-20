function [ReachableSet] = GetReachableSetAt( varargin )
	%Description:
	%	Gets the reachable set at time t if the POB_Feedback object contains enough
	%	information.
	%
	%Usage:
	%	[ ReachableSet ] = pob1.GetReachableSetAt( t , word_index )
	%	[ ReachableSet ] = pob1.GetReachableSetAt( t , word_index , 'PwT' , PwT )
	
	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	[ pob1 , t , word_index , Px0 , L , PwT ] = InputProcessing_GetReachableSetAt( varargin{:} );

	lcsas0 = pob1.System;
	[ n_x , n_u , n_y , n_w , n_v ] = lcsas0.Dimensions();

	% Create MPC Matrices w.r.t. the desired word.
	target_word = L.words{word_index};
	TimeHorizon = length(target_word);
	[S_w_i,S_u_i,~,J_i,f_bar_i] = get_mpc_matrices(lcsas0,'word',target_word);

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	% Find a Valid Knowledge Sequence In The Controller Such That Contains target_word at time t.
	num_knowl_sequences = size(pob1.PossibleLSequences,2);
	matching_seq_id = -1;
	for knowl_seq_index = 1:num_knowl_sequences
		temp_lang_at_t = pob1.PossibleLSequences(t,knowl_seq_index);
		if temp_lang_at_t.contains(target_word)
			matching_seq_id = knowl_seq_index;
		end
	end
	if matching_seq_id == -1
		error(['There were no matches for the desired language at time ' num2str(t) '.' ])
	end

	% Use That Gain
	K = pob1.F_set{matching_seq_id}; k = pob1.u0_set{matching_seq_id};
	K_trimmed = K([1:n_u*t],[1:n_w*TimeHorizon]);
	k_trimmed = k([1:n_u*t],1);


	S_w_trimmed = S_w_i([1:n_x*(t+1)],[1:n_w*TimeHorizon]);
	S_u_trimmed = S_u_i([1:n_x*(t+1)],[1:n_u*t]);
	J_trimmed = J_i([1:n_x*(t+1)],:);
	f_bar_prime = S_w_trimmed*f_bar_i;

	x = ( S_w_trimmed + S_u_trimmed * K_trimmed)*PwT + S_u_trimmed*k_trimmed + J_trimmed * Px0 + f_bar_prime ;

	ReachableSet = x.projection(n_x*t+[1:n_x]);

end

function [ pob_out , t , word_index , X0 , L , PwT ] = InputProcessing_GetReachableSetAt( varargin )
	%Description:
	%	Require the system to use the version of LCSAS with X0 defined.

	%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Default Assignments %%
	%%%%%%%%%%%%%%%%%%%%%%%%%

	pob_out = varargin{1};
	t = varargin{2};
	word_index = varargin{3};

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Identifying any extra inputs %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if nargin > 3
		argidx = 4;
		while argidx <= nargin
			switch varargin{argidx}
				case 'PwT'
					PwT = varargin{argidx+1};
					argidx = argidx + 2;
				otherwise
					error(['Unexpected input to GetReachableSetAt(): ' varargin{argidx} ])
			end
		end
	end

	%%%%%%%%%%%%%%%%%%%%
	%% Input Checking %%
	%%%%%%%%%%%%%%%%%%%%

	%% Check that the POB_Feedback object has a valid X0 object inside of it.

	if isempty(pob_out.System.X0)
		error('The input system does not have a well-defined X0 field. Please define it first.')
	end

	allowed_X0_types = {'double','Polyhedron'};
	X0_type_comparison = strcmp(class(pob_out.System.X0),allowed_X0_types);
	if ~any( X0_type_comparison )
		error(['X0 is not of the right type. It must be either a double or Polyhedron.'])
	end

	X0 = pob_out.System.X0;

	%% Language

	L = pob_out.System.L;

	%% Check that the mode index is allowable according to the language L

	if (1 > word_index) || (word_index > L.cardinality())
		error(['There is not a word ' num2str(word_index) ' in the language L. It has only ' num2str(L.cardinality()) ' words.'  ])
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Create Any Extra Variables %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% PwT: A polyhedron containing the disturbance sequence w = [ w(1); w(2); ... ]
	if ~exist('PwT')
		PwT = 1;
		for mode_val = L.words{word_index}(1:t)
			PwT = PwT * pob_out.System.Dyn(mode_val).P_w;
		end
	end

end 