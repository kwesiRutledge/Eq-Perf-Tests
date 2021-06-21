function [ H_cl , h_cl ] = get_closed_loop_consistent_internal_behavior_set_matrices( lcsas_in , H_ol , h_ol , x0 , K , k , path_in )
	%Description
	%	Receives as input: 
	%		(i) the open loop consistent internal behavior matrices, 
	%		(ii) a desired controller parameterized by K and k
	%	and returns the matrices H and h which correspond to the closed loop consistent internal behavior set (which is a polytope).
	%
	%Assumptions:
	%	Note that it is assumed that (H_ol,h_ol) is given in the dimension n_x*(T+1)+(n_u+n_w)*T.

	%% Input Processing %%

	if ~isvector(path_in)
		error(['The input path_in must be a vector. Instead, it''s size is ' num2str(size(path_in)) '.'])
    end

    [ lcsas_in , path_in , TimeHorizon , magic_t ] = input_processing_gclcibsm( lcsas_in , H_ol , h_ol , x0 , K , k , path_in );
    
	%% Constants %%

	[ n_x , n_u , n_y , n_w , n_v ] = lcsas_in.Dimensions();

	[ S_w , S_u , ~ , J , f_bar ] = lcsas_in.get_mpc_matrices('All Words');

	last_L = path_in(end);
	lL_card = last_L.cardinality();
    
    % magic_t = TimeHorizon - 1; %The last time acknowledgeable by the algorithm    

	%% Algorithm %%

	nonK_prefactor = {}; K_prefactor = {};
	h_independent_factor = {}; h_dependent_factor = {};

	H_cl = []; h_cl = [];

	for tll_index = 1:last_L.cardinality()
 		[~,word_id] = lcsas_in.L.contains(last_L.words{tll_index});

 		nonK_prefactor{tll_index} = ...
 			[ [ zeros(n_x*(magic_t+1),(tll_index-1)*n_w*magic_t), S_w{word_id}, zeros(n_x*(magic_t+1),(lL_card-tll_index)*n_w*magic_t) ] ;
				zeros(n_u*magic_t,lL_card*n_w*magic_t) ;
				eye(lL_card*n_w*magic_t) ;
				zeros(lL_card*n_x,lL_card*n_w*magic_t) ];

		K_prefactor{tll_index} = [	S_u{word_id} ;
							eye(n_u*magic_t);
							zeros(lL_card*n_w*magic_t,n_u*magic_t);
							zeros(lL_card*n_x,n_u*magic_t) ] ;

		%This matrix is needed to select the "w" trajectory corresponding to the current mode.
		SelectWMatrix = [ zeros(n_w*magic_t,(tll_index-1)*n_w*magic_t), eye(n_w*magic_t), zeros(n_w*magic_t,(lL_card-tll_index)*n_w*magic_t) ];

		H_cl = [	H_cl ;
					H_ol * ( nonK_prefactor{tll_index} + K_prefactor{tll_index} * K([1:n_u*magic_t],[1:n_w*magic_t]) * SelectWMatrix ) ] ;
		
		h_independent_factor{tll_index} = h_ol - ...
			H_ol * [ J{word_id}*x0 + S_w{word_id}*f_bar{word_id} ;
					zeros(n_u*magic_t,1) ;
					zeros(lL_card*n_w*magic_t,1) ; 
					repmat(x0,lL_card,1) ];

		h_dependent_factor{tll_index} = - ...
			H_ol * [ S_u{word_id} ;
					eye(n_u*magic_t) ;
					zeros(lL_card*n_w*magic_t,n_u*magic_t) ; 
					zeros(lL_card*n_x,n_u*magic_t) ];

		h_cl = [ 	h_cl ;
					h_independent_factor{tll_index} + h_dependent_factor{tll_index}*k([1:n_u*magic_t],1) ];

		% % Create Constraints
		% dummy_var_bound_constraints = dummy_var_bound_constraints + [ -temp_K_bound <= temp_dummy_eta{end} <= temp_K_bound ] + [ temp_dummy_y{end} <= temp_K_bound ];

	end

end

function [ lcsas_in , path_in , TimeHorizon , magic_t ] = input_processing_gclcibsm( lcsas_in , H_ol , h_ol , x0 , K , k , path_in ) 
    %Description:
    %
    
    %% Constants
    L = lcsas_in.L;
    
    %% 
    
    %% Trim Language
    trim_length = length(L.words{1}) + 1 - length(path_in);

    L_trimmed = L.trim_by(trim_length);
    
    %% Replace Path
    new_path = [];
    for path_index = 1:length(path_in)
        %Create new path_element;
        new_path = [ new_path ; path_in(path_index).trim_by(trim_length) ];
    end
    
    
    %% Create Outputs
    
    lcsas_in.L = L_trimmed;
    path_in = new_path;
    TimeHorizon = length(L.words{1});
    magic_t = length(path_in) - 1;
    
end