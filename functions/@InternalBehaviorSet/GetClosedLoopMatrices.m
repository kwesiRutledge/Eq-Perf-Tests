function [ A_cl , b_cl , Ae_cl , be_cl ] = GetClosedLoopMatrices( ibs )
	%Description
	%	Receives as input: 
	%		(i) the open loop consistent internal behavior matrices, 
	%		(ii) a desired controller parameterized by K and k
	%	and returns the matrices H and h which correspond to the closed loop consistent internal behavior set (which is a polytope).
	%
	%Assumptions:
	%	Note that it is assumed that (H_ol,h_ol) is given in the dimension n_x*(T+1)+(n_u+n_w)*T.
	%
	%Notes:
	%	We will assume that there is only one gain applied and that it is for disturbances coming from the FIRST
	%	word in the last language.
	%	This file was mostly copied from get_closed_loop_consistent_internal_behavior_set_matrices.m

    % [ lcsas_in , path_in , TimeHorizon , t ] = input_processing_gclcibsm( lcsas_in , H_ol , h_ol , x0 , K , k , path_in );

    %% Input Processing %%
    ibs_settings = ibs.ibs_settings;
    if ~strcmp(ibs_settings.fb_type,'state')
    	error('The InternalBehaviorSet must have feedback type ''state'' in order to be computed effectively.')
    end

	%% Constants %%

	System = ibs.System;
	t = ibs.t;
	KnowledgeSequence = ibs.KnowledgeSequence;

	[ n_x , n_u , n_y , n_w , n_v ] = System.Dimensions();
	[ S_w , S_u , ~ , J , f_bar ] = System.get_mpc_matrices('All Words');

	last_L = KnowledgeSequence(end);
	lL_card = last_L.cardinality();
    
	A_ol = ibs.A; b_ol = ibs.b;
	Ae_ol = ibs.Ae; be_ol = ibs.be;

	K = ibs.K; k = ibs.k;

	%This matrix is needed to select the "w" trajectory corresponding to the current mode.
	SelectWMatrix1 = ibs.SelectW();

	%% Algorithm %%

	nonK_prefactor = []; K_prefactor = [];
	h_independent_factor = []; h_dependent_factor = [];

	H_cl = []; h_cl = [];

	Ae_p = [zeros(n_u*t,n_x*(t+1)),-eye(n_u*t),zeros(n_u*t,ibs.Dim-n_x*(t+1)-(n_u)*t)] + ...
			K([1:n_u*t],[1:n_w*t]) * SelectWMatrix1;
	be_p = [-k([1:n_u*t],1)];

	%% Create Outputs %%

	A_cl = A_ol; b_cl = b_ol;
	Ae_cl = [ Ae_ol ; Ae_p ];
	be_cl = [ be_ol ; be_p ];

end

function [ lcsas_in , path_in , TimeHorizon , t ] = input_processing_gclcibsm( lcsas_in , H_ol , h_ol , x0 , K , k , path_in ) 
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
    t = length(path_in) - 1;
    
end