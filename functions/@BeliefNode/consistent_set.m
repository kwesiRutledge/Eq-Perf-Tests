function [Consist_set, full_set ] = consistent_set(varargin)
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
	%
	%Inputs:
	%	lcsas 	- An LCSAS object. Hopefully the dimensions are all appropriately checked so that
	%			  the state is the proper size in all Aff_Dyn(). Inputs as well, etc.
	%	t 		- The time of interest
	%	L 		- The set of words under consideration.
	%			  We would like to find the set of states for which it is possible to reach when under ALL switching
	%			  sequences in this set with the same inputs (albeit with different disturbances).
	%
	%Example Usage:
	%	[Consist_set, full_set] = BN.consistency_set(lcsas,P_u,P_x0,'fb_method','output')
	%	[Consist_set, full_set] = BN.consistency_set(lcsas,P_u,P_x0,'fb_method','state')
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
	BN = varargin{1};
	lcsas = varargin{2};
	P_u = varargin{3};
	P_x0 = varargin{4};
	
	if ~isa(lcsas,'LCSAS')
		error('Expecting the first input to be a LCSAS object.')
	end

	varargin_idx = 5;
	while varargin_idx <= nargin
		switch varargin{varargin_idx}
			case 'fb_method'
				fb_type = varargin{varargin_idx+1};
				if ~(strcmp(fb_type,'state') || strcmp(fb_type,'output'))
					error(['Invalid feedback type: ' fb_type ])
				end
				varargin_idx = varargin_idx + 2;
			otherwise
				error('Unexpected additional input.')
		end
	end

	if ~exist('fb_type')
		fb_type = 'state';
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

	L = lcsas.L;

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
	I_blockx = []; I_blockx2 = []; I_blocky = [];
	C_block = []; Cv_block = [];
	for word_ind = 1:length(L.words)
		H_block(end+[1:size(Hc{word_ind},1)],end+[1:size(Bwc{word_ind},2)]) = -Hc{word_ind}*Bwc{word_ind};
		S_block(end+[1:size(Sc{word_ind},1)],[1:size(Sc{word_ind},2)]) = -Sc{word_ind};
		J_block(end+[1:size(Jc{word_ind},1)],end+[1:size(Jc{word_ind},2)]) = -Jc{word_ind};
		f_block(end+[1:size(Hc{word_ind}*fc{word_ind},1)],1) = Hc{word_ind}*fc{word_ind};

		I_blockx(end+[1:n_x*(t+1)],[1:n_x*(t+1)]) = eye(n_x*(t+1));
		I_blocky(end+[1:n_y*(t+1)],[1:n_y*(t+1)]) = eye(n_y*(t+1));

		C_block(end+[1:n_y*(t+1)],end+[1:n_x*(t+1)]) = -[Cc{word_ind} ; zeros(n_y,n_x*t), lcsas.Dyn( L.words{word_ind}(t+1) ).C ];
		Cv_block(end+[1:n_y*(t+1)],end+[1:n_v*(t+1)]) = -[Cvc{word_ind},zeros(size(Cvc{word_ind},1),n_v);zeros(n_y,size(Cvc{word_ind},2)), lcsas.Dyn( L.words{word_ind}(t+1) ).C_v ];
		I_blockx2(end+[1:n_x*(t+1)],end+[1:n_x*(t+1)]) = eye(n_x*(t+1));
	end

	%% Constructing the Sets
	if strcmp(fb_type,'state')
	    
		P_eta = P_uT * P_wT * P_x0_L;

		% Create the Equality Constraints
		Hc = {}; Sc = {}; Jc = {}; fc = {};
		for word_ind = 1:length(L.words)
			[Hc{word_ind},Sc{word_ind},~,Jc{word_ind},fc{word_ind}] = lcsas.get_mpc_matrices('word',L.words{word_ind}(1:t));
		end

		%Create the set of feasible (x,u,w,x0) tuples
		full_set = Polyhedron('A',[zeros(size(P_eta.A,1),n_x*(t+1)),P_eta.A],'b',P_eta.b,'Ae',[I_blockx, S_block, H_block, J_block],'be',f_block );

		%Project the above set to create the set of feasible observed trajectories (x,u)
		Consist_set = [ eye(n_x*(t+1) + n_u*t), zeros(n_x*(t+1) + n_u*t, length(L.words)*(n_w*t + n_x) ) ] * full_set;
	
	else
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
    							'Ae',[zeros(size(S_block,1),size(I_blocky,2)),S_block, H_block, zeros(size(S_block,1),size(Cv_block,2)), J_block, I_blockx2; ...
    								  I_blocky, zeros(size(I_blocky,1),size(S_block,2)+size(H_block,2)), Cv_block , zeros(size(I_blocky,1),size(J_block,2)) , C_block ], ...
    							'be', [f_block;zeros(size(I_blocky,1),1)] );

		%Project the above set to create the set of feasible observed trajectories (x,u)
		Consist_set = [ eye(n_y*(t+1) + n_u*t), zeros(n_x*(t+1) + n_u*t, length(L.words)*(n_w*t + n_v*(t+1) + n_x + n_x*(t+1)) ) ] * full_set;

	end
		

end