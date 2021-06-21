function varargout = internal_behavior_set( varargin )
	%internal_behavior_set.m
	%Description:
	%	Finds a polyhedron that describes what pairs of state, input, and disturbance sequences are compatible/feasible from ALL
	%	switching sequences defined by the sequence KnowledgeSequence.
	%
	%	The KnowledgeSequence is a vector describing the following sequence of estimates about the true state of the system;
	%		KnowledgeSequence = [ L_{0} , L_{1} , ... , L_{T-1} ]
	%
	%	To place the system into clearer focus. We have a Language-Constrained Switched Affine System (LCSAS):
	%
	%	x_{t+1} = A_{q_t} x_t + B_{q_t} u_t + f_{q_t} + w_t
	%	
	%	where:
	%			- q_t is a natural number that describes the current mode at time t
	%			- w_t belongs to the set W_{q_t} which varies with the mode
	%	The consistency set can also be written as
	%
	%									  { [x]  | }
	%		IB_lcsas(KnowledgeSequence) = { [u]  | }
	%									  { [w]  | }
	%									  { [x0] | }
	%		
	%					{ [y]  | }
	%					{ [u]  | }
	%		C(\sigma) = { [w]  | }
	%					{ [v]  | }
	%					{ [x0] | }
	%					{ [x]  | }
	%
	%Inputs:
	%	lcsas 		- An array of Aff_Dyn() objects. Hopefully the dimensions are all appropriately checked so that
	%				  the state is the proper size in all Aff_Dyn(). Inputs as well, etc.
	%	use_proj 	- Boolean (true or false).
	%				  Used to tell the function to either skip the creation of Consist_set (false) or complete the
	%				  computation of Consist_set (true) which requires projection operations to be called and may be very slow.
	%
	%Example Usage:
	%	[ internal_behavior_set ] = lcsas.internal_behavior_set( KnowledgeSequence )
	%	[ internal_behavior_set ] = lcsas.internal_behavior_set( KnowledgeSequence , 'fb_type' , 'state' )
	%
	%Assumptions:
	%	This formulation assumes that the system does not include a disturbed measurements. i.e. We can perfectly observe the state
	%	at each time.

	%% Input Processing
	[ lcsas0 , KnowledgeSequence , ibs_settings ] = ip_internal_behavior_set( varargin{:} );

	%% Constants

	KS_len = length(KnowledgeSequence);
	T = length(lcsas0.L.words{1});
	% L_Tm1 = KnowledgeSequence(end);
    
    P_u = lcsas0.U;

	%% Algorithm
        
    % Create Internal Behavior Set Using Each of the ClairvoyantSequences
    ibs = [];
        
    ib_components_at_ = {};
    % ib_components_at_{1} = lcsas0.X0;

    for t = 1:KS_len %KnowledgeSequence estimates L_0, ... , L_{T-2}
        L_t = KnowledgeSequence(t); 

        %Collect the ib_at_t sets for each word in L_t
        ib_components_at_{t} = {};
        for word_idx = 1:L_t.cardinality()
            L_t_projected = Language(L_t.words{word_idx});
            temp_ib_sets(word_idx) = internal_behavior_set_at_t( lcsas0 , t , L_t_projected , ibs_settings );

            % Transform the set into the dimension described by max_card.
            if word_idx == 1
                ib_components_at_{t} = modify_ibs_with_respect_to_T( lcsas0 , KnowledgeSequence , temp_ib_sets(word_idx) , L_t_projected , t , T );
            else
                ib_components_at_{t} = ib_components_at_{t}.intersect( ...
                    modify_ibs_with_respect_to_T( lcsas0 , KnowledgeSequence , temp_ib_sets(word_idx) , L_t_projected , t , T ) ...
                );	
            end
        end
    end

    ib_at_T = ib_components_at_{1};
    for t = 2:length(ib_components_at_)
        ib_at_T = ib_at_T.intersect(ib_components_at_{t});
    end
    
    ibs = ib_at_T ;
    
	%%%%%%%%%%%%%%%%%%%%
	%% Create Outputs %%
	%%%%%%%%%%%%%%%%%%%%

	varargout{1} = ibs;


end

function [ lcsas0 , KnowledgeSequence , ibs_settings ] = ip_internal_behavior_set( varargin )
	%Description:
	%	Input processing function for internal_behavior_set().
	%	Verifies that there:
	%	- are at least four inputs,
	%	- LCSAS input checks,

	%% Check Argument Numbers
	expected_number_of_arguments = 2;

	if nargin < expected_number_of_arguments
		error(['Expect at least four arguments.'])
	end

	%% First Required Input: lcsas
	lcsas0 = varargin{1};
	lcsas0.check('X0','U'); %Verify that X0 and U exist and are both Polyhedron objects.

	%% Second Required Input: KnowledgeSequence
	KnowledgeSequence = varargin{2};
	if ~isa(KnowledgeSequence,'Language')
		error(['KnowledgeSequence  is not of type ''Language''. Instead it is of type ' class(KnowledgeSequence) '.'])
	end

	if ~isvector(KnowledgeSequence)
		error(['KnowledgeSequence must be a vector or scalar of Language objects.'])
	end

	%% Define Defaults
	ibs_settings = struct('fb_type','state','reduce_flag',false);

	all_fb_types = {'state','output'};

	%% Check the additional inputs
	argument_index = expected_number_of_arguments+1;
	while argument_index <= nargin
		switch varargin{argument_index}
			case 'fb_type'
				ibs_settings.fb_type = varargin{argument_index+1};
				if ~any(strcmp(ibs_settings.fb_type,all_fb_types))
					error(['Unexpected fb_type value: ' ibs_settings.fb_type ])
				end
				argument_index = argument_index + 2;
			case 'reduce_flag'
				ibs_settings.reduce_flag = varargin{argument_index+1};
				argument_index = argument_index + 2;
			otherwise
				error(['Unexpected input to internal_behavior_set: ' varargin{argument_index} ])
		end
	end
end

function [ ib_at_t ] = internal_behavior_set_at_t( lcsas0 , t , L_t , ibs_settings )
	%Description:
	%	Creates the set of behaviors that are consistent with the hypothesis L_t at time t.
	%	This set does not consider hypotheses at previous times and is thus a sometimes over approximation
	%	of what can be expected.
	%
	%Notes:
	%	It's possible that this can be simplified. We are enforcing a lot of constraints that don't need to be enforced.
	%	i.e. the set is enforcing constraints also on times t-1,t-2, etc. but it does not need to.
	%
	%	Each trajectory contained in the internal behavior set starts from t=0. So, the internal behavior_set_at_t=1
	%	will contain 
	%
	%			[ x0 ]
	%		x = [ x1 ] , u = [ u0 ] , w = [ w0 ]
	%

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	[ n_x , n_u , n_y , n_w , n_v ] = lcsas0.Dimensions();

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	% Create large disturbance set from Cartesian products
	P_wT = 1; P_x0_L = 1;
	for word_idx = 1:L_t.cardinality()
		for symb_idx = 1:t
			P_wT = P_wT * lcsas0.Dyn( L_t.words{word_idx}(symb_idx) ).P_w;
		end
		%Initial Condition Set for Each Word in L
		P_x0_L = P_x0_L * lcsas0.X0;
    end
    
    % Created Disturbance Sets
    P_uT = 1;
    for t_idx = 1:t
        P_uT = P_uT * lcsas0.U;
    end

    % Create mpc matrices for each word in the language L
	Hc = {}; Sc = {}; Jc = {}; fc = {}; Cc = {}; Bwc = {}; Cvc = {};
	for word_ind = 1:length(L_t.words)
		[Hc{word_ind},Sc{word_ind},Cc{word_ind},Jc{word_ind},fc{word_ind},Bwc{word_ind},Cvc{word_ind}] = lcsas0.get_mpc_matrices('word',L_t.words{word_ind}(1:t));
	end

	H_block = []; S_block = []; J_block = []; f_block = [];
	I_blockx = []; I_blockx2 = [];  I_blocky = [];
	C_block = []; Cv_block = [];
	for word_ind = 1:length(L_t.words)
		H_block(end+[1:size(Hc{word_ind},1)],end+[1:size(Bwc{word_ind},2)]) = Hc{word_ind}*Bwc{word_ind};
		S_block(end+[1:size(Sc{word_ind},1)],[1:size(Sc{word_ind},2)]) = Sc{word_ind};
		J_block(end+[1:size(Jc{word_ind},1)],end+[1:size(Jc{word_ind},2)]) = Jc{word_ind};
		f_block(end+[1:size(Hc{word_ind}*fc{word_ind},1)],1) = Hc{word_ind}*fc{word_ind};

		I_blockx(end+[1:n_x*(t+1)],[1:n_x*(t+1)]) = eye(n_x*(t+1));
		I_blocky(end+[1:n_y*(t+1)],[1:n_y*(t+1)]) = eye(n_y*(t+1));

	end

	%% Constructing the Sets
	if strcmp(ibs_settings.fb_type,'state')
	    
		P_eta = P_uT * P_wT * P_x0_L;

		%Create the set of feasible (x,u,w,x0) tuples
		ib_at_t = Polyhedron(	'A',[zeros(size(P_eta.A,1),n_x*(t+1)),P_eta.A],'b',P_eta.b, ...
								'Ae',[-I_blockx, S_block, H_block, J_block],'be',-f_block );

		% Add Equality Constraints Between Words in Internal Behavior Sets
		for word_ind = 2:L_t.cardinality()
			Ltc = L_t.cardinality();
			matching_x_block = [ eye(n_x*(t+1)) , zeros( n_x*(t+1) , (-2+word_ind)*n_x*(t+1) ) , -eye(n_x*(t+1)) , zeros(n_x*(t+1) , (Ltc-word_ind)*n_x*(t+1) ) ];
			ib_at_t.Ae = [ 	ib_at_t.Ae ;
							matching_x_block , zeros(n_x*(t+1),size(S_block,2)) , zeros(n_x*(t+1),size(H_block,2)) , zeros(n_x*(t+1),size(J_block))  ];

			ib_at_t.be = [ ib_at_t.be ; zeros(n_x*(t+1),1) ];
		end

	else strcmp(ibs_settings.fb_type,'output')

		%Also introduce the measurement disturbance into the equation
		P_vT = 1; 
		for word_idx = 1:length(L.words)
			for symb_idx = 1:t+1
				P_vT = P_vT * lcsas.Dyn( L.words{word_idx}(symb_idx) ).P_v;
			end
    	end

    	% Introduce Sets
    	for word_ind = 1:L.cardinality()
    		C_block(end+[1:n_y*(t+1)],end+[1:n_x*(t+1)]) = [Cc{word_ind} ; zeros(n_y,n_x*t), lcsas.Dyn( L.words{word_ind}(t+1) ).C ];
			Cv_block(end+[1:n_y*(t+1)],end+[1:n_v*(t+1)]) = [Cvc{word_ind},zeros(size(Cvc{word_ind},1),n_v);zeros(n_y,size(Cvc{word_ind},2)), lcsas.Dyn( L.words{word_ind}(t+1) ).C_v ];
			I_blockx2(end+[1:n_x*(t+1)],end+[1:n_x*(t+1)]) = eye(n_x*(t+1));
		end

    	P_eta = P_uT * P_wT * P_vT * P_x0_L;

    	%Create the set of feasible (x,u,w,x0) tuples
    	ib_at_t = Polyhedron(	'A',[zeros(size(P_eta.A,1),n_y*(t+1)),P_eta.A,zeros(size(P_eta.A,1),length(L.words)*n_x*(t+1))],'b',P_eta.b, ...
    							'Ae',[zeros(size(S_block,1),size(I_blocky,2)),S_block, H_block, zeros(size(S_block,1),size(Cv_block,2)), J_block, -I_blockx2; ...
    								  I_blocky, zeros(size(I_blocky,1),size(S_block,2)+size(H_block,2)), -Cv_block , zeros(size(I_blocky,1),size(J_block,2)) , -C_block ], ...
    							'be', [-f_block;zeros(size(I_blocky,1),1)] );
    end

    if (~ib_at_t.isEmptySet) && ibs_settings.reduce_flag
		ib_at_t.minHRep; %Make sure to do this to simplify some of the future projections.
	end

end

function [ ibs_extended ] = modify_ibs_with_respect_to_T( lcsas0 , KnowledgeSequence , ibs_in , L_in , t , T )
	%Description:
	%	Modifies the internal behavior set defined for a time t and map it into the dimension
	%	of the expected behaviors at time T.
	%

	%% Input Checking

	% Verify that the input language has a single word in it.
	if L_in.cardinality() ~= 1
		error(['modify_ibs_with_respect_to_T() expects a language with only one word inside of it.'])
	end

	%% Constants

	[ n_x , n_u , n_y , n_w , n_v ] = lcsas0.Dimensions();

    % trimmedKnowledgeSequence = KnowledgeSequence([2:end]);
    
	[ max_card , max_card_index ] = KnowledgeSequence.find_maximum_cardinality_in_sequence();

	ibs_e_dim = n_x * (T+1) + n_u*T + n_w*T*max_card + n_x*max_card;

	%% Algorithm

	%Find the location of the input Language L_in relative to the language at max_card_t
	L_mct = KnowledgeSequence(max_card_index);
	[tf,L_in_index] = L_mct.contains(L_in.words{1});
 
	if ~tf
		error(['The input Language L_in for modify_ibs_with_respect_to_T() is not contained within L_mct.'])
	end

	L_in_index_as_binary = zeros(1,L_mct.cardinality());
	L_in_index_as_binary( L_in_index ) = 1;

	A_Prefactor = [ ...
		eye( n_x * (t+1) ), zeros( n_x*(t+1) , ibs_e_dim - n_x * (t+1) ) ;
		zeros( n_u*(t) , n_x*(T+1) ), eye(n_u*t) , zeros( n_u*(t) , ibs_e_dim - n_x*(T+1) - n_u*t ) ;
		zeros( n_w*(t) , n_x*(T+1) + n_u*T ) , kron( L_in_index_as_binary , [ eye(n_w*t) , zeros(n_w*t,n_w*(T-t)) ] ) , zeros(n_w*t, ibs_e_dim - n_x*(T+1) - n_u*T - n_w*T*max_card) ;
		zeros( n_x , ibs_e_dim - n_x*max_card), kron( L_in_index_as_binary , eye(n_x) ) ...
		];

	ibs_extended = Polyhedron( ...
		'A' , ibs_in.A *A_Prefactor, 'b' , ibs_in.b , ...
		'Ae', ibs_in.Ae*A_Prefactor, 'be', ibs_in.be );

end 