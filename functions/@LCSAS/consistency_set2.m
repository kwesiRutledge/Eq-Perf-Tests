function [Consist_set, full_set ] = consistency_set2(varargin)
	%consistent_set.m
	%Description:
	%	Finds a polyhedron that describes what pairs of states and input sequences are compatible/feasible from ALL
	%	switching sequences defined by L.
	%
	%	To place the system into clearer focus. We have a Language-Constrained Switched Affine System (LCSAS):
	%
	%	x_{t+1} = A_{q_t} x_t + B_{q_t} u_t + f_{q_t} + w_t
	%	
	%	where:
	%			- q_t is a natural number that describes the current mode at time t
	%			- w_t belongs to the set W_{q_t} which varies with the mode
	%	The consistency set can also be written as"
	%
	%							{ [y]  | 			[y]  }
	%		full_set(\sigma) =	{ [u]  | \exists 	[u]  }
	%							{  					[w]  }
	%							{  					[x0]  }
	%
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
	%	[Phi_t_L] = consistent_set(lcsas,t,L)
	%	[Consist_set, full_set] = consistent_set(lcsas,t,L,P_u,P_x0)
	%	[Consist_set, full_set] = consistent_set(lcsas,t,L,P_u,P_x0,'fb_method','state')
	%	[Consist_set, full_set] = consistent_set(lcsas,t,L,P_u,P_x0,'fb_method','state','use_proj',false)
	%
	%Assumptions:
	%	This formulation assumes that the system does not include a disturbed measurements. i.e. We can perfectly observe the state
	%	at each time.

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

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

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	n_x = size(lcsas.Dyn(1).A,2);
	n_u = size(lcsas.Dyn(1).B,2);
	n_w = size(lcsas.Dyn(1).B_w,2);
	n_y = size(lcsas.Dyn(1).C,1);
	n_v = size(lcsas.Dyn(1).C_v,2);

	if ~exist('fb_type')
		fb_type = 'state';
	end

	if ~exist('use_proj')
		use_proj = true;
	end

	if ~exist('reduce_flag')
		reduce_flag = true;
	end

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	%Create large disturbance set from Cartesian products
	P_wT = 1; P_x0_L = 1;
	for word_idx = 1:length(L.words)
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

		C_block(end+[1:n_y*(t+1)],end+[1:n_x*(t+1)]) = [Cc{word_ind} ; zeros(n_y,n_x*t), lcsas.Dyn( L.words{word_ind}(t+1) ).C ];
		%Cv_block(end+[1:n_y*(t+1)],end+[1:n_v*(t+1)]) = [Cvc{word_ind},zeros(size(Cvc{word_ind},1),n_v);zeros(n_y,size(Cvc{word_ind},2)), lcsas.Dyn( L.words{word_ind}(t+1) ).C_v ];
		I_blockx2(end+[1:n_x*(t+1)],end+[1:n_x*(t+1)]) = eye(n_x*(t+1));
	end

	% I_block_x0 = []; 
	% block_select_x0 = zeros(L.cardinality()*n_x,size(I_blockx2));
	% for word_ind = 1:L.cardinality()
	% 	block_select_x0(end+[1:n_x],:) = [ zeros(n_x,(n_x*(t+1))*(word_ind-1)) eye(n_x)  ]
	% end

	%% Constructing the Sets
	if strcmp(fb_type,'state')
	    
		error('Feedback type is unsupported at the moment.')

		% P_eta = P_uT * P_wT * P_x0_L;

		% %Create the set of feasible (x,u,w,x0) tuples
		% full_set = Polyhedron(	'A',[zeros(size(P_eta.A,1),n_x*(t+1)),P_eta.A],'b',P_eta.b, ...
		% 						'Ae',[-I_blockx, S_block, H_block, J_block],'be',-f_block );

	else
		%Also introduce the measurement disturbance into the equation
		P_vT = 1; 
		for word_idx = 1:length(L.words)
			for symb_idx = 1:t+1
				P_vT = P_vT * lcsas.Dyn( L.words{word_idx}(symb_idx) ).P_v;
			end
    	end

    	P_eta = P_uT * P_wT * P_x0_L;

    	%Create the set of feasible (x,u,w,x0) tuples
    	full_set = Polyhedron(	'A',[	zeros(size(P_eta.A,1),n_y*(t+1)) , P_eta.A;
    									P_vT.A*[ I_blocky , -C_block*S_block , -C_block*H_block,-C_block*J_block ] ], ...
    							'b',[ 	P_eta.b ; P_vT.b + P_vT.A*C_block*f_block ]	);

	end

	if (~full_set.isEmptySet) && reduce_flag
		full_set.minHRep; %Make sure to do this to simplify some of the future projections.
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Constructing the Consistency Sets %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if use_proj

		if strcmp(fb_type,'state')

			%Project the above set to create the set of feasible observed trajectories (x,u)
			Consist_set = [ eye(n_x*(t+1) + n_u*t), zeros(n_x*(t+1) + n_u*t, length(L.words)*(n_w*t + n_x) ) ] * full_set;
		
		else

	    	if ~full_set.isEmptySet
				%Project the above set to create the set of feasible observed trajectories (x,u)
				%Consist_set = [ eye(n_y*(t+1) + n_u*t), zeros(n_x*(t+1) + n_u*t, length(L.words)*(n_w*t + n_v*(t+1) + n_x + n_x*(t+1)) ) ] * full_set;
				%Consist_set = full_set.affineMap([ eye(n_y*(t+1) + n_u*t), zeros(n_x*(t+1) + n_u*t, length(L.words)*(n_w*t + n_v*(t+1) + n_x + n_x*(t+1)) ) ],'vrep')
				Consist_set = full_set.projection([1:n_y*(t+1) + n_u*t]);
			else
				Consist_set = Polyhedron('A',[ [1;-1], zeros(2,n_y*(t+1) + n_u*t-1) ],'b',[1;-2]);
			end
		end
		%Consist_set.minHRep; %Used to make sure that future projections are simpler to compute.

	else
		Consist_set = Polyhedron('A',[[1;-1],zeros(2,n_x*(t+1) + n_u*t-1)], ...
								 'b',[1;-2]);
	end
		

end