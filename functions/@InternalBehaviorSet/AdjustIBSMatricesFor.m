function [ A , b , Ae , be ] = adjust_ibs_matrices_for( ibs , t0 , L0 ,  A0 , b0 , Ae0, be0 )
	%Description:
	%	Modifies the internal behavior set defined for a time t and map it into the dimension
	%	of the expected behaviors at time T.
	%

	%% Input Checking

	% Verify that the input language has a single word in it.
	if L0.cardinality() ~= 1
		error(['modify_ibs_with_respect_to_T() expects a language with only one word inside of it.'])
	end

	%% Constants

	KnowledgeSequence = ibs.KnowledgeSequence;
	System = ibs.System;

	[ n_x , n_u , n_y , n_w , n_v ] = System.Dimensions();

	T = length(System.L.words{1});
	t = t0;

    % trimmedKnowledgeSequence = KnowledgeSequence([2:end]);
    
	[ max_card , max_card_index ] = KnowledgeSequence.find_maximum_cardinality_in_sequence();

	ibs_e_dim = n_x * (T+1) + n_u*T + n_w*T*max_card + n_x*max_card;

	%% Algorithm

	%Find the location of the input Language L0 relative to the language at max_card_t
	L_mct = KnowledgeSequence(max_card_index);
	[tf,L0_index] = L_mct.contains(L0.words{1});
 
	if ~tf
		error(['The input Language L0 for modify_ibs_with_respect_to_T() is not contained within L_mct.'])
	end

	L0_index_as_binary = zeros(1,L_mct.cardinality());
	L0_index_as_binary( L0_index ) = 1;

	A_Prefactor = [ ...
		eye( n_x * (t+1) ), zeros( n_x*(t+1) , ibs_e_dim - n_x * (t+1) ) ;
		zeros( n_u*(t) , n_x*(T+1) ), eye(n_u*t) , zeros( n_u*(t) , ibs_e_dim - n_x*(T+1) - n_u*t ) ;
		zeros( n_w*(t) , n_x*(T+1) + n_u*T ) , kron( L0_index_as_binary , [ eye(n_w*t) , zeros(n_w*t,n_w*(T-t)) ] ) , zeros(n_w*t, ibs_e_dim - n_x*(T+1) - n_u*T - n_w*T*max_card) ;
		zeros( n_x , ibs_e_dim - n_x*max_card), kron( L0_index_as_binary , eye(n_x) ) ...
		];

	A = A0 * A_Prefactor;
	b = b0;

	Ae = Ae0 * A_Prefactor;
	be = be0;

	% ibs_extended = Polyhedron( ...
	% 	'A' , ibs_in.A *A_Prefactor, 'b' , ibs_in.b , ...
	% 	'Ae', ibs_in.Ae*A_Prefactor, 'be', ibs_in.be );

end 