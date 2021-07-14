function [ A , b , Ae , be ] = AdjustIBSMatricesFor( ibs , tau , L0 ,  A0 , b0 , Ae0, be0 )
	%Description:
	%	Modifies the internal behavior set defined for a time t and map it into the dimension
	%	of the expected behaviors at time T.
	%

	%% Input Checking

	%% Constants

	KnowledgeSequence = ibs.KnowledgeSequence;
	System = ibs.System;
	ibs_settings = ibs.ibs_settings;

	[ n_x , n_u , n_y , n_w , n_v ] = System.Dimensions();

	t = ibs.t;

    % trimmedKnowledgeSequence = KnowledgeSequence([2:end]);
    
	[ max_card , max_card_index ] = KnowledgeSequence.find_maximum_cardinality_in_sequence();

	%% Algorithm

	%Find the location of the input Language L0 relative to the language at max_card_t
	L_mct = KnowledgeSequence(max_card_index);

	L0_indicies = -ones(L_mct.cardinality(),1);
	L0_indicies_as_binary = zeros(L_mct.cardinality(),1);

	mat_L0 = [];
	for word_index = 1:L0.cardinality()
		[tf,L0_index] = L_mct.contains(L0.words{word_index});
		if tf
			temp_row = zeros(L_mct.cardinality(),1)';
			temp_row(L0_index) = 1;

			mat_L0 = [ mat_L0 ; temp_row ];

			L0_indicies(L0_index) = 1;

			% L0_indicies_as_binary(L0_index) = 1;
		end
	end
 
	if ~tf
		error(['The input Language L0 for modify_ibs_with_respect_to_T() is not contained within L_mct.'])
	end

	% L0_index_as_binary = zeros(1,L_mct.cardinality());
	% L0_index_as_binary( L0_index ) = 1;

	%mat_L0 = diag(L0_indicies_as_binary);

	switch ibs_settings.fb_type
	case 'state'
		ibs_e_dim = n_x * (t+1) + n_u*t + n_w*t*max_card + n_x; %Expected Dimension of Expected Behavior Set

		A_Prefactor = [ ...
			eye( n_x * (tau+1) ), zeros( n_x*(tau+1) , ibs_e_dim - n_x * (tau+1) ) ;
			zeros( n_u*(tau) , n_x*(t+1) ), eye(n_u*tau) , zeros( n_u*(tau) , ibs_e_dim - n_x*(t+1) - n_u*tau ) ;
			zeros( n_w*(tau)*L0.cardinality() , n_x*(t+1) + n_u*t ) , kron( mat_L0 , [ eye(n_w*tau) , zeros(n_w*tau,n_w*(t-tau)) ] ) , zeros(n_w*tau*L0.cardinality(), ibs_e_dim - n_x*(t+1) - n_u*t - n_w*t*max_card) ;
			zeros( n_x , ibs_e_dim - n_x), eye(n_x) ...
			];

	case 'output'
		ibs_e_dim = n_y * (t+1) + n_u*t + n_w*t*max_card + (n_v+n_x)*(t+1)*max_card + n_x*max_card; %Expected Dimension of Expected Behavior Set

		A_Prefactor = [ ...
			eye( n_y * (tau+1) ), zeros( n_y*(tau+1) , ibs_e_dim - n_y * (tau+1) ) ;
			zeros( n_u*(tau) , n_y*(t+1) ), eye(n_u*tau) , zeros( n_u*(tau) , ibs_e_dim - n_y*(t+1) - n_u*tau ) ;
			zeros( n_w*(tau)*L0.cardinality() , n_y*(t+1) + n_u*t ) , kron( mat_L0 , [ eye(n_w*tau) , zeros(n_w*tau,n_w*(t-tau)) ] ) , zeros(n_w*tau*L0.cardinality(), ibs_e_dim - n_y*(t+1) - n_u*t - n_w*t*max_card) ;
			zeros( n_v*(tau+1)*L0.cardinality() , n_y*(t+1) + n_u*t + n_w*t*max_card ) , kron( mat_L0 , [ eye(n_v*(tau+1)) , zeros(n_v*(tau+1),n_v*(t-tau)) ] ) , zeros( n_v*(tau+1)*L0.cardinality() , ibs_e_dim - n_y*(t+1) - n_u*t - n_w*t*max_card - n_v*(t+1)*max_card ) ;
			zeros( n_x*L0.cardinality() , n_y*(t+1) + n_u*t + (n_w*t+n_v*(t+1))*max_card), kron( mat_L0 , eye(n_x) ), zeros( n_x*L0.cardinality() , ibs_e_dim - n_y*(t+1) - n_u*t - n_w*t*max_card - n_v*(t+1)*max_card - n_x*max_card );
			zeros( n_x*(tau+1)*L0.cardinality() , ibs_e_dim - n_x*(t+1)*L0.cardinality() ) , kron( mat_L0 , [ eye(n_x*(tau+1)) , zeros(n_x*(tau+1),n_x*(t-tau)) ] )
			];

	otherwise
		error(['Unexpected fb_type value in AdjustIBSMatricesFor().'])
	end

	A = A0 * A_Prefactor;
	b = b0;

	Ae = Ae0 * A_Prefactor;
	be = be0;

	% ibs_extended = Polyhedron( ...
	% 	'A' , ibs_in.A *A_Prefactor, 'b' , ibs_in.b , ...
	% 	'Ae', ibs_in.Ae*A_Prefactor, 'be', ibs_in.be );

end 