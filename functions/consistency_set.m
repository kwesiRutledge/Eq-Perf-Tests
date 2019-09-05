function [Phi_t_L] = consistency_set(varargin)
%consistency_set.m
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
%	lcsas 	- An array of Aff_Dyn() objects. Hopefully the dimensions are all appropriately checked so that
%			  the state is the proper size in all Aff_Dyn(). Inputs as well, etc.
%	t 		- The time of interest
%	L 		- The set of words under consideration.
%			  We would like to find the set of states for which it is possible to reach when under ALL switching
%			  sequences in this set with the same inputs (albeit with different disturbances).
%
%Example Usage:
%	[Phi_t_L] = consistency_set(lcsas,t,L)
%	[Phi_t_L] = consistency_set(lcsas,t,L,P_u,P_x0)
%
%Assumptions:
%	This formulation assumes that the system does not include a disturbed measurements. i.e. We can perfectly observe the state
%	at each time.

%% Input Processing

switch nargin
	case 3
		error('3 Argument version coming soon?')
	case 5
		lcsas = varargin{1};
		t = varargin{2};
		L = varargin{3};
		P_u = varargin{4};
		P_x0 = varargin{5};
	otherwise
		error('Check the number of input arguments. Too many or too few given.')
end


if (t < 0)
	error(['t must have a value greater than 0.'])
end

for word_ind = 1:length(L)
	if t > length(L{word_ind})
		error('t should not be larger than any of the words in L.')
	end
end

n_x = size(lcsas(1).A,2);

%% Algorithm
height.P_w_blocks = 0;
for word_ind = 1:length(L)
	Pw_A_blocks{word_ind} = []; Pw_b_blocks{word_ind} = [];
	for sys_ind = L{word_ind}(1:t)
		%Create each block of the A matrices
		Pw_A_blocks{word_ind}(end+[1:size(lcsas(sys_ind).P_w.A,1)],end+[1:size(lcsas(sys_ind).P_w.A,2)]) = lcsas(sys_ind).P_w.A;
		Pw_b_blocks{word_ind}(end+[1:size(lcsas(sys_ind).P_w.A,1)],1) = lcsas(sys_ind).P_w.b;
	end
	%Now create the full A and b matrices
	height.P_w_blocks = height.P_w_blocks + size(lcsas(sys_ind).P_w.A,1);
end

temp_A = [zeros(height.P_w_blocks,n_x)];
Pw_A_big_block = []; Pw_b_big_block = [];
Pu_A_big_block = []; Pu_b_big_block = [];
Px0_A_big_block = []; Px0_b_big_block = [];

for word_ind = 1:length(L)
	Pw_A_big_block(end+[1:size(Pw_A_blocks{word_ind},1)],end+[1:size(Pw_A_blocks{word_ind},2)]) = Pw_A_blocks{word_ind};
	Pw_b_big_block(end+[1:size(Pw_b_blocks{word_ind},1)],1) = Pw_b_blocks{word_ind};

	Px0_A_big_block(end+[1:size(P_x0.A,1)],end+[1:n_x]) = P_x0.A;
	Px0_b_big_block(end+[1:size(P_x0.A,1)],1) = P_x0.b;
end

Pu_A_big_block = kron(eye(t),P_u.A);
Pu_b_big_block = kron(ones(t,1),P_u.b);

height = [];
[height.Pw_A_bb, width0.Pw_A_bb] = size(Pw_A_big_block);
[height.Pu_A_bb, width0.Pu_A_bb] = size(Pu_A_big_block);
[height.Px0_A_bb, width0.Px0_A_bb] = size(Px0_A_big_block);

temp_A = [	zeros(height.Pu_A_bb,n_x*(t+1)), 								Pu_A_big_block,	zeros(height.Pu_A_bb,width0.Pw_A_bb+width0.Px0_A_bb);
			zeros(height.Pw_A_bb,n_x*(t+1)+width0.Pu_A_bb),					Pw_A_big_block,	zeros(height.Pw_A_bb,width0.Px0_A_bb);
			zeros(height.Px0_A_bb,n_x*(t+1)+width0.Pu_A_bb+width0.Pw_A_bb),	Px0_A_big_block];

temp_b = [ Pu_b_big_block; Pw_b_big_block; Px0_b_big_block ];

% Create the Equality Constraints
rt = [zeros(n_x,n_x*t) eye(n_x)];

Hc = {}; Sc = {}; Jc = {}; fc = {};
for word_ind = 1:length(L)
	[Hc{word_ind},Sc{word_ind},~,Jc{word_ind},fc{word_ind}] = get_mpc_matrices(lcsas,'word',L{word_ind}(1:t));
end

H_block = []; S_block = []; J_block = []; f_block = []; I_block = [];
for word_ind = 1:length(L)
	H_block(end+[1:size(Hc{word_ind},1)],end+[1:size(Hc{word_ind},2)]) = -Hc{word_ind};
	S_block(end+[1:size(Sc{word_ind},1)],[1:size(Sc{word_ind},2)]) = -Sc{word_ind};
	J_block(end+[1:size(Jc{word_ind},1)],end+[1:size(Jc{word_ind},2)]) = -Jc{word_ind};
	f_block(end+[1:size(Hc{word_ind}*fc{word_ind},1)],1) = Hc{word_ind}*fc{word_ind};

	I_block(end+[1:n_x*(t+1)],[1:n_x*(t+1)]) = eye(n_x*(t+1));
end;

%Create test_poly2
Phi_t_L = Polyhedron('A',temp_A,'b',temp_b,'Ae',[I_block, S_block, H_block, J_block],'be',f_block );

end