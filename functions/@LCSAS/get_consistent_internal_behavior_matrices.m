function [ H_Phi, h_Phi ] = get_consistent_internal_behavior_matrices(varargin)
	%get_consistent_internal_behavior_matrices.m
	%Description:
	%	Defines the matrices which define the internal behavior set Phi containing the measurements, inputs, disturbances,
	%	initial conditions, and state trajectories that are consistent with a given mode language L for the system lcsas.
	%
	%	To place the system into clearer focus. We have a Language-Constrained Switched Affine System (LCSAS):
	%
	%	x_{t+1} = A_{q_t} x_t + B_{q_t} u_t + f_{q_t} + w_t
	%	
	%	where:
	%			- q_t is a natural number that describes the current mode at time t
	%			- w_t belongs to the set W_{q_t} which varies with the mode
	%	The internal behavior set can also be written as"
	%						{ [y]  |  y_t = C_{q_t}x_t + v_t 	 }
	%						{ [u]  |  						 	 }
	%		Phi(\sigma) =	{ [w]  |  w_t \in lcsas.Dyn(q_t).P_w }
	%						{ [v]  |  v_t \in lcsas.Dyn(q_t).P_v }
	%						{ [x0] |  x_0 \in P_x0 }
	%						{ [x]  |  x_t = A_{q_t} x_t + v_t 	 }
	%	and q_t comes from a sequence in L.
	%
	%Inputs:
	%	lcsas 		- An array of Aff_Dyn() objects. Hopefully the dimensions are all appropriately checked so that
	%				  the state is the proper size in all Aff_Dyn(). Inputs as well, etc.
	%	t 			- The time of interest
	%	L 			- The set of words under consideration.
	%				  We would like to find the set of states for which it is possible to reach when under ALL switching
	%				  sequences in this set with the same inputs (albeit with different disturbances).
	%	use_proj 	- Boolean (true or false).
	%				  Used to tell the function to either skip the creation of Consist_set (false) or complete the
	%				  computation of Consist_set (true) which requires projection operations to be called and may be very slow.
	%
	%Example Usage:
	%	[ H_Phi , h_Phi ] = lcsas.get_consistent_internal_behavior_matrices( t , L , x0 , P_u , P_x0 )
	%
	%Assumptions:
	%	This formulation assumes that the system does not include a disturbed measurements. i.e. We can perfectly observe the state
	%	at each time.

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	[ lcsas , t , L , P_u , P_x0 , fb_type , use_proj , reduce_flag ] = ip_get_consistent_internal_behavior_matrices(varargin{:});

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	[ n_x , n_u , n_y , n_w , n_v ] = lcsas.Dimensions();

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	%Create large disturbance set from Cartesian products
	P_wT = 1; P_x0_L = 1;
	for word_idx = 1:L.cardinality()
		for symb_idx = 1:t
			P_wT = P_wT * lcsas.Dyn( L.words{word_idx}(symb_idx) ).P_w;
		end
		%Initial Condition Set for Each Word in L
		P_x0_L = P_x0_L * P_x0;
    end
    %Created Disturbance Sets

    P_uT = 1;
    for t_idx = 1:t
        P_uT = P_uT * P_u;
    end

    % Create mpc matrices for each word in the language L
	Hc = {}; Sc = {}; Jc = {}; fc = {}; Cc = {}; Bwc = {}; Cvc = {};
	for word_ind = 1:length(L.words)
		[Hc{word_ind},Sc{word_ind},Cc{word_ind},Jc{word_ind},fc{word_ind},Bwc{word_ind},Cvc{word_ind}] = lcsas.get_mpc_matrices('word',L.words{word_ind}(1:t));
	end

	H_block = []; S_block = []; J_block = []; f_block = [];
	I_blockx = []; I_blockx2 = [];  I_blocky = [];
	C_block = []; Cv_block = [];
	for word_ind = 1:length(L.words)
		H_block(end+[1:size(Hc{word_ind},1)],end+[1:size(Bwc{word_ind},2)]) = Hc{word_ind}*Bwc{word_ind};
		S_block(end+[1:size(Sc{word_ind},1)],[1:size(Sc{word_ind},2)]) = Sc{word_ind};
		J_block(end+[1:size(Jc{word_ind},1)],end+[1:size(Jc{word_ind},2)]) = Jc{word_ind};
		f_block(end+[1:size(Hc{word_ind}*fc{word_ind},1)],1) = Hc{word_ind}*fc{word_ind};

		I_blockx(end+[1:n_x*(t+1)],[1:n_x*(t+1)]) = eye(n_x*(t+1));
		I_blocky(end+[1:n_y*(t+1)],[1:n_y*(t+1)]) = eye(n_y*(t+1));

	end

	% I_block_x0 = []; 
	% block_select_x0 = zeros(L.cardinality()*n_x,size(I_blockx2));
	% for word_ind = 1:L.cardinality()
	% 	block_select_x0(end+[1:n_x],:) = [ zeros(n_x,(n_x*(t+1))*(word_ind-1)) eye(n_x)  ]
	% end

	%% Constructing the Matrices
	if strcmp(fb_type,'state')
	    
		P_eta = P_uT * P_wT * P_x0_L;

		%Create the set of feasible (x,u,w,x0) tuples
		H_Phi = [zeros(size(P_eta.A,1),n_x*(t+1)),P_eta.A];
		H_Phi = [ 	H_Phi ;
					[-I_blockx, S_block, H_block, J_block] ;
					-[-I_blockx, S_block, H_block, J_block] ];

		h_Phi = P_eta.b;
		h_Phi = [ 	h_Phi ;
					-f_block ;
					-(-f_block) ]; %Adding equality constraints

	else

		error('This function is not ready yet.')

		for word_ind = 1:length(L.words)
			C_block(end+[1:n_y*(t+1)],end+[1:n_x*(t+1)]) = [Cc{word_ind} ; zeros(n_y,n_x*t), lcsas.Dyn( L.words{word_ind}(t+1) ).C ];
			Cv_block(end+[1:n_y*(t+1)],end+[1:n_v*(t+1)]) = [Cvc{word_ind},zeros(size(Cvc{word_ind},1),n_v);zeros(n_y,size(Cvc{word_ind},2)), lcsas.Dyn( L.words{word_ind}(t+1) ).C_v ];
			I_blockx2(end+[1:n_x*(t+1)],end+[1:n_x*(t+1)]) = eye(n_x*(t+1));
		end

		%Also introduce the measurement disturbance into the equation
		P_vT = 1; 
		for word_idx = 1:length(L.words)
			for symb_idx = 1:t+1
				P_vT = P_vT * lcsas.Dyn( L.words{word_idx}(symb_idx) ).P_v;
			end
    	end

    	P_eta = P_uT * P_wT * P_vT * P_x0_L;

    	%Create the set of feasible (x,u,w,x0) tuples
    	full_set = Polyhedron(	'A',[zeros(size(P_eta.A,1),n_y*(t+1)),P_eta.A,zeros(size(P_eta.A,1),length(L.words)*n_x*(t+1))],'b',P_eta.b, ...
    							'Ae',[zeros(size(S_block,1),size(I_blocky,2)),S_block, H_block, zeros(size(S_block,1),size(Cv_block,2)), J_block, -I_blockx2; ...
    								  I_blocky, zeros(size(I_blocky,1),size(S_block,2)+size(H_block,2)), -Cv_block , zeros(size(I_blocky,1),size(J_block,2)) , -C_block ], ...
    							'be', [-f_block;zeros(size(I_blocky,1),1)] );

	end

end

function [ lcsas , t , L , P_u , P_x0 , fb_type , use_proj , reduce_flag ] = ip_get_consistent_internal_behavior_matrices(varargin)
	%Description:
	%	Processes the input given to get_consistent_internal_behavior_matrices() according to the
	%	expectations listed in the overall function.

	%%%%%%%%%%%%%%%%%%%%%
	%% Handling Inputs %%
	%%%%%%%%%%%%%%%%%%%%%

	if nargin < 3
		error('Not enough input arguments.')
	end

	lcsas = varargin{1};
	t = varargin{2};
	L = varargin{3};
	P_u = varargin{4};
	P_x0 = varargin{5};
	
	if ~isa(lcsas,'LCSAS')
		error('Expecting the first input to be a LCSAS object.')
	end

	if ~isa(L,'Language')
		error('Expecting the language input to be a Language object.')
	end

	varargin_idx = 6;
	while varargin_idx <= nargin
		switch varargin{varargin_idx}
			case 'fb_method'
				fb_type = varargin{varargin_idx+1};
				if ~(strcmp(fb_type,'state') || strcmp(fb_type,'output'))
					error(['Invalid feedback type: ' fb_type ])
				end
				varargin_idx = varargin_idx + 2;
			case 'use_proj'
				use_proj = varargin{varargin_idx+1};
				if ~islogical( use_proj )
					error('The flag for ''use_proj'' should be a boolean.')
				end
				varargin_idx = varargin_idx + 2;
			case 'reduce_representation'
				reduce_flag = varargin{varargin_idx+1};
				if ~islogical( reduce_flag )
					error('The flag for ''reduce_flag'' should be a boolean.')
				end
				varargin_idx = varargin_idx + 2;
			otherwise
				error('Unexpected additional input.')
		end
	end

	if (t < 0)
		error(['t must have a value greater than 0.'])
	end

	for word_ind = 1:length(L.words)
		if t > length(L.words{word_ind})
			error('t should not be larger than any of the words in L.')
		end
	end

	%%%%%%%%%%%%%%%%%%
	%% Set Defaults %%
	%%%%%%%%%%%%%%%%%%

	if ~exist('fb_type')
		fb_type = 'state';
	end

	if ~exist('use_proj')
		use_proj = true;
	end

	if ~exist('reduce_flag')
		reduce_flag = true;
	end

end