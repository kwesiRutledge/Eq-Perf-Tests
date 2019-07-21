%experiment42.m
%% Description
%	The purpose of this experiment is to:
%	- test some of the new modifications to the Aff_Dyn() object.

%% Initialization
% include_fcns

clear all; close all; clc;

%% Constants
h = 0.1; %Discretization step
epsil0 = 0.6; %Coefficient of restitution
m = 1;

%Define the Piecewise Affine System
dim = 2;

A1 = [ 1, h; 0 1 ];
B1 = [0; (h/m)];
C1 = eye(dim);
f1 = zeros(dim,1);

eta_v = 0; eta_w = 0.2;
Pv1 = Polyhedron('lb',-eta_v*ones(1,dim) ,'ub',eta_v*ones(1,dim));
Pw1 = Polyhedron('lb',-eta_w*ones(1,dim) ,'ub',eta_w*ones(1,dim));

ad1 = Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1);

A2 = [1,0;0,-epsil0];

ad2 = Aff_Dyn(A2,B1,f1,C1,Pw1,Pv1);

hs1 = [ad1,ad2];

clear ad1 ad2

%Define a simpler Hybrid System
A1 = eye(dim);
B1 = eye(dim);
C1 = eye(dim);
f1 = [0;1];

ad1 = Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1);

f2 = [1;0];
f3 = -f1;
f4 = -f2;

hs2 = [	ad1,...
		Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1),...
		Aff_Dyn(A1,B1,f3,C1,Pw1,Pv1),...
		Aff_Dyn(A1,B1,f4,C1,Pw1,Pv1)];

% Select Matrix
RT = @(n,t,T) [zeros(n,n*t) eye(n) zeros(n,n*(T-t))];

%% Testing the mpc matrix generator on Switched System 1
[H1,S1,~,J1,~] = get_mpc_matrices(hs1,'word',[1,2]);
[H,S,~,f_bar] = get_mpc_matrices(hs1,'word',[1,1,1,1,2,2,2,2]);

P_x0 = Polyhedron('lb',-ones(1,dim),'ub',ones(1,dim));
P_u  = Polyhedron('lb',-2,'ub',2); 

X = H1*(Pw1*Pw1) + S1*(P_u*P_u) + J1*P_x0;

figure;
plot([zeros(dim,2*dim),eye(dim)]*X,'color','magenta')

%% Testing the Visualization of Switched System 2
clear H1 S1 J1
clear Hc Sc Jc

clear P_u P_x0

%Define Sets
eta_u = 0; eta_x0 = 0.2;
P_u = Polyhedron('lb',-eta_u*ones(1,dim) ,'ub',eta_u*ones(1,dim));
P_x0 = Polyhedron('lb',-eta_x0*ones(1,dim),'ub',eta_x0*ones(1,dim));

Hc = {}; Sc = {}; Jc = {}; fc = {};

[Hc{1},Sc{1},~,Jc{1},fc{1}] = get_mpc_matrices(hs2,'word',[1,1,2]);
[Hc{2},Sc{2},~,Jc{2},fc{2}] = get_mpc_matrices(hs2,'word',[1,2,2]);
[Hc{3},Sc{3},~,Jc{3},fc{3}] = get_mpc_matrices(hs2,'word',[1,2,4]);

T = 3;

figure;
for seq_ind = 1:length(Hc)
	X = (Hc{seq_ind}*(Pw1*Pw1*Pw1+fc{seq_ind})) + Sc{seq_ind}*(P_u*P_u*P_u) + Jc{seq_ind}*P_x0;
	subplot(1,length(Hc),seq_ind)
	hold on;
	plot(RT(dim,0,3)*X,'Color','white')
	plot(RT(dim,1,3)*X,'Color','cyan')
	plot(RT(dim,2,3)*X,'Color','magenta')
	plot(RT(dim,3,3)*X,'Color','orange')
	axis([-5 5 -5 5])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Testing the Definition of the 'Detected Sets' %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear H1 S1 J1
clear Hc Sc Jc

clear P_u P_x0

%Define Sets
eta_u = 0.25; eta_x0 = 0.5;
P_u = Polyhedron('lb',-eta_u*ones(1,dim) ,'ub',eta_u*ones(1,dim));
P_x0 = Polyhedron('lb',-eta_x0*ones(1,dim),'ub',eta_x0*ones(1,dim));

Hc = {}; Sc = {}; Jc = {}; fc = {};

[Hc{1},Sc{1},~,Jc{1},fc{1}] = get_mpc_matrices(hs2,'word',[1,1,2]);
[Hc{2},Sc{2},~,Jc{2},fc{2}] = get_mpc_matrices(hs2,'word',[1,2,2]);
[Hc{3},Sc{3},~,Jc{3},fc{3}] = get_mpc_matrices(hs2,'word',[1,2,4]);

figure;
for seq_ind = 1:length(Hc)
	X = (Hc{seq_ind}*(Pw1*Pw1*Pw1+fc{seq_ind})) + Sc{seq_ind}*(P_u*P_u*P_u) + Jc{seq_ind}*P_x0;
	subplot(1,length(Hc),seq_ind)
	hold on;
	plot(RT(dim,0,3)*X,'Color','white')
	plot(RT(dim,1,3)*X,'Color','cyan')
	plot(RT(dim,2,3)*X,'Color','magenta')
	plot(RT(dim,3,3)*X,'Color','orange')
	axis([-5 5 -5 5])
end

%There should now be inputs that create similar final states:
test_poly = Polyhedron( 'A',[ 	zeros(T*(size(Pw1.A,1)),dim) kron(eye(T),Pw1.A) zeros(T*(size(Pw1.A,1)),2*T*dim+2*dim) ; ...
								zeros(T*(size(Pw1.A,1)),dim+T*dim) kron(eye(T),Pw1.A) zeros(T*(size(Pw1.A,1)),T*dim+2*dim) ; ...
								zeros(T*(size(P_u.A,1)),dim+2*T*dim) kron(eye(T),P_u.A) zeros(T*(size(P_u.A,1)),2*dim) ; ...
								zeros((size(P_x0.A,1)),dim+3*T*dim) P_x0.A zeros((size(P_x0.A,1)),dim); ...
								zeros((size(P_x0.A,1)),dim+3*T*dim+dim) P_x0.A ], ...
						'b',[	kron(ones(2*T,1),Pw1.b); ...
								kron(ones(T,1),P_u.b); ...
								kron(ones(2,1),P_x0.b)], ...
						'Ae',[	eye(dim),	-RT(dim,3,3)*Hc{2},				zeros(size(-RT(dim,3,3)*Hc{3})),-RT(dim,3,3)*Sc{2},-RT(dim,3,3)*Jc{2},		zeros(size(-RT(dim,3,3)*Jc{3})); ...
								eye(dim),	zeros(size(-RT(dim,3,3)*Hc{2})),-RT(dim,3,3)*Hc{3},				-RT(dim,3,3)*Sc{3},zeros(size(-RT(dim,3,3)*Jc{2})),-RT(dim,3,3)*Jc{3}], ...
						'be', [RT(dim,3,3)*Hc{2}*fc{2}; ...
								RT(dim,3,3)*Hc{3}*fc{3}]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Implementing a Function for generating the Belief Tree %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear n

%First we'll need a function that generates the test_poly function above.

%Inputs to such a function would be:

L1 = {[1,1,2]}; % A Language
t1 = 3;			% The desired time we want to observe.
hs_ut = hs2;			%The hybrid system

%Algorithm
n1 = size(hs1(1).A,2);
height.P_w_blocks = 0;

for word_ind = 1:length(L1)
	Pw_A_blocks{word_ind} = []; Pw_b_blocks{word_ind} = [];
	for sys_ind = L1{word_ind}(1:t1)
		%Create each block of the A matrices
		Pw_A_blocks{word_ind}(end+[1:size(hs_ut(sys_ind).P_w.A,1)],end+[1:size(hs_ut(sys_ind).P_w.A,2)]) = hs_ut(sys_ind).P_w.A;
		Pw_b_blocks{word_ind}(end+[1:size(hs_ut(sys_ind).P_w.A,1)],1) = hs_ut(sys_ind).P_w.b;
	end
	%Now create the full A and b matrices
	height.P_w_blocks = height.P_w_blocks + size(hs_ut(sys_ind).P_w.A,1);
end

temp_A = [zeros(height.P_w_blocks,n1)];
Pw_A_big_block = []; Pw_b_big_block = [];
Pu_A_big_block = []; Pu_b_big_block = [];
Px0_A_big_block = []; Px0_b_big_block = [];

for word_ind = 1:length(L1)
	Pw_A_big_block(end+[1:size(Pw_A_blocks{word_ind},1)],end+[1:size(Pw_A_blocks{word_ind},2)]) = Pw_A_blocks{word_ind};
	Pw_b_big_block(end+[1:size(Pw_b_blocks{word_ind},1)],end+[1:size(Pw_b_blocks{word_ind},2)]) = Pw_b_blocks{word_ind};

	Px0_A_big_block(end+[1:size(P_x0.A,1)],end+[1:n1]) = P_x0.A;
	Px0_b_big_block(end+[1:size(P_x0.A,1)],1) = P_x0.b;
end

Pu_A_big_block = kron(eye(t1),P_u.A);
Pu_b_big_block = kron(ones(t1,1),P_u.b);

[height.Pw_A_bb, width0.Pw_A_bb] = size(Pw_A_big_block);
[height.Pu_A_bb, width0.Pu_A_bb] = size(Pu_A_big_block);
[height.Px0_A_bb, width0.Px0_A_bb] = size(Px0_A_big_block);


temp_A = [	zeros(height.Pw_A_bb,n1*(t1+1)), 								Pu_A_big_block,	zeros(height.Pu_A_bb,width0.Pw_A_bb+width0.Px0_A_bb);
			zeros(height.Pu_A_bb,n1*(t1+1)+width0.Pu_A_bb),					Pw_A_big_block,	zeros(height.Pw_A_bb,width0.Px0_A_bb);
			zeros(height.Px0_A_bb,n1*(t1+1)+width0.Pw_A_bb+width0.Pu_A_bb),	Px0_A_big_block];

temp_b = [ Pu_b_big_block; Pw_b_big_block; Px0_b_big_block ];

% Create the Equality Constraints
rt = [zeros(n1,n1*t1) eye(n1)];

Hc = {}; Sc = {}; Jc = {}; fc = {};
for word_ind = 1:length(L1)
	[Hc{word_ind},Sc{word_ind},~,Jc{word_ind},fc{word_ind}] = get_mpc_matrices(hs_ut,'word',L1{word_ind}(1:t1));
end

H_block = []; S_block = []; J_block = []; f_block = []; I_block = [];
for word_ind = 1:length(L1)
	H_block(end+[1:size(Hc{word_ind},1)],end+[1:size(Hc{word_ind},2)]) = -Hc{word_ind};
	S_block(end+[1:size(Sc{word_ind},1)],[1:size(Sc{word_ind},2)]) = -Sc{word_ind};
	J_block(end+[1:size(Jc{word_ind},1)],end+[1:size(Jc{word_ind},2)]) = -Jc{word_ind};
	f_block(end+[1:size(Hc{word_ind}*fc{word_ind},1)],1) = Hc{word_ind}*fc{word_ind};

	I_block(end+[1:n1*(t1+1)],end+[1:n1*(t1+1)]) = eye(n1*(t1+1));
end;

%Create test_poly2
test_poly2 = Polyhedron('A',temp_A,'b',temp_b,'Ae',[I_block, S_block, H_block, J_block],'be',f_block );

%Create test_poly2 using the function
test_poly3 = consistency_set(hs2,t1,L1,P_u,P_x0);

disp('Belief Tree Components: Consistency Set (Pt. 1)')
disp(' ')
disp('Are test_poly2 and test_poly3 the same polytope?')
disp(['test_poly2 contains test_poly3? ' num2str(test_poly2.contains(test_poly3))])
disp(['test_poly3 contains test_poly2? ' num2str(test_poly3.contains(test_poly2))])
disp(' ')

%Define the Predecessor Operator?
t2 = 1;
n_u = size(P_u.A,2);
Phi1 = consistency_set(hs2,t2,{[1,2,2],[1,2,4]},P_u,P_x0);
Phi2 = consistency_set(hs2,t2+1,{[1,2,2],[1,2,4]},P_u,P_x0);

Projxu_Phi1 = [eye((t2+1)*n1+t2*n_u) zeros((t2+1)*n1+t2*n_u,size(Phi1.A,2)-((t2+1)*n1+t2*n_u))]*Phi1;
Projxu_Phi2 = [[eye((t2+1)*n1);zeros(t2*n_u,(t2+1)*n1)], zeros((t2+1)*n1+t2*n_u,n1),[zeros((t2+1)*n1,t2*n_u);eye(t2*n_u)],  zeros((t2+1)*n1+t2*n_u,size(Phi2.A,2)-((t2+1+1)*n1+t2*n_u))]*Phi2;

temp_intersx = Projxu_Phi1.intersect(Projxu_Phi2);

disp('Belief Tree Components: Intersections of Consistency Sets (Pt. 2)')
disp(' ')
disp(['Is Projxu_Phi1 \cap Projxu_Phi2 empty? ' num2str(temp_intersx.isEmptySet)] )
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%
%% Belief Tree Test %%
%%%%%%%%%%%%%%%%%%%%%%

disp('Belief Tree Assembly (Test 1)')
disp(' ')

%Define inputs
L2 = {[1,1,2],[1,2,2],[1,2,4]};
T2 = 3;
hs2;

%Constants
n2 = size(hs2(1).A,1);
m2 = size(hs2(1).B,2);

%Create first node
node0.subset = L2;
node0.t = 0;
% node0.parent = NaN;
node1.subset = [1:length(L2)];
node1.t = 0;

c_level = [node0];

Phi0 = P_x0*P_u;

nodes0 = [node0];
nodes1 = [node1]
edges0 = [];
edges1 = []; %Numerical version of the edges matrix

for tau = 1:T2
	%Each belief will be indexed by a time. (i.e. I hold X belieft at time t)
	for node_ind = 1:length(c_level) %Iterate through all nodes that are stored in the c_level array
		%Get All Combinations of the node's subset
		node_p_set = {};
		for comb_length = 1:length(c_level(node_ind).subset)
			temp_combs = nchoosek([1:length(c_level(node_ind).subset)],comb_length);
			for comb_ind = 1:size(temp_combs,1)
				node_p_set{end+1} = temp_combs(comb_ind,:);
			end
		end

		%Build Phi Sets
		Phi_sets = {}; visible_transitions = [];
		for p_set_ind = 1:length(node_p_set)
			Phi_sets{p_set_ind} = consistency_set(hs2,tau,{c_level(node_ind).subset{node_p_set{p_set_ind}}},P_u,P_x0);
			%If any of the Phi's are empty,
			%then it is impossible for a transition to exist between the node c_level(node_ind) and the node associated with Phi
			if ~Phi_sets{p_set_ind}.isEmptySet
				visible_transitions = [visible_transitions p_set_ind];
			end
		end 

		%For each possible transition, see if its transition set is completely contained by another transition set
		for ind_ut = 1:length(node_p_set)
			Projx_Phi1 = [eye(n2*(tau+1)+m2*tau) zeros(n2*(tau+1)+m2*tau,Phi_sets{ind_ut}.Dim-n2*(tau+1)-m2*tau)]*Phi_sets{ind_ut};
			%vt_tilde = visible_transitions(visible_transitions ~= ind_ut);
			for ind_ch = [ind_ut+1:length(node_p_set)]
				disp(['ind_ch = ' num2str(ind_ch)])
				Projx_Phi2 = [eye(n2*(tau+1)+m2*tau) zeros(n2*(tau+1)+m2*tau,Phi_sets{ind_ch}.Dim-n2*(tau+1)-m2*tau)]*Phi_sets{ind_ch};
				temp_diff = Projx_Phi1 \ Projx_Phi2;
				if (temp_diff.isEmptySet)
					disp(['temp_diff.isEmptySet = ' num2str(temp_diff.isEmptySet) ' for:'])
					disp(['- Phi1(' num2str([c_level(node_ind).subset{node_p_set{ind_ut}}]) ')' ])
					disp(['- Phi2(' num2str([c_level(node_ind).subset{node_p_set{ind_ch}}]) ')' ])
					disp(' ')
					visible_transitions(ind_ut) = ind_ch;
				% elseif (ind_ch == ind_ut) && (~Projx_Phi1.isEmptySet)
				% 	visible_transitions(ind_ut) = ind_ch;
				end 
			end
		end
		visible_transitions = unique(visible_transitions);
		
		% Create edges and the next level of the tree
		for edge_ind = 1:length(visible_transitions)
			temp_node.subset = {c_level(node_ind).subset{node_p_set{visible_transitions(edge_ind)}}}; 
			temp_node.t = tau;
			
			temp_node1.subset = [];
			for subset_ind = 1:length(temp_node.subset)
				for L2_ind = 1:length(L2)
					if all(L2{L2_ind} == temp_node.subset{subset_ind})
						temp_node1.subset = [temp_node1.subset, L2_ind];
					end
				end
			end
			temp_node1.t = tau;
			%Add temporary node to the nodes list
			node_in_set_already = false;
			for clevel_ind = 1:length(nodes0)
				% disp(['language_eq(nodes0(clevel_ind).subset,temp_node.subset) = ' num2str(language_eq(nodes0(clevel_ind).subset,temp_node.subset)) ])
				if (language_eq(nodes0(clevel_ind).subset,temp_node.subset)) && (nodes0(clevel_ind).t == temp_node.t)
					node_in_set_already = true;
				end
			end
			if ~node_in_set_already
				nodes0 = [nodes0,temp_node];
				nodes1 = [nodes1,temp_node1];
			end
			%Create edge using this new node and add it to the edges list
			temp_edge = [c_level(node_ind),temp_node];
			edges0 = [edges0;temp_edge];
		end

	end

	%Create next level of the tree
	c_level = [];
	for node_ind = 1:length(nodes0)
		if nodes0(node_ind).t == tau
			c_level = [c_level,nodes0(node_ind)];
		end
	end

	disp(['There are ' num2str(length(c_level)) ' nodes at time tau = ' num2str(tau) '.' ])
end

%Clean Up Tree

for vert_ind = 1:length(nodes0)
	temp_node = nodes0(vert_ind);
	for edge_ind = 1:size(edges0,1)
		%Place a vertex number instead of a whole vertex into the node location in the edge matrix for a more compact representation
		if language_eq(edges0(edge_ind,1).subset,temp_node.subset) && (edges0(edge_ind,1).t == temp_node.t)
			edges1(edge_ind,1) = vert_ind;
		elseif language_eq(edges0(edge_ind,2).subset,temp_node.subset) && (edges0(edge_ind,2).t == temp_node.t)
			edges1(edge_ind,2) = vert_ind;
		end
	end
end

%% Enumerate All Paths through the Belief Tree

bTree.V = nodes0;
bTree.E = edges0;

disp(' ')
disp(['Find the different belief sequences that are possible. A belief language.'])
disp(' ')

%Start from the top.
b_seq = {};
b_seq{1} = {L2};
for seq_time = 2:length(L2{1})
	for pref_ind = 1:length(b_seq{seq_time-1})
		temp_prefix = b_seq{seq_time-1}{pref_ind};
		%Add suffixes based on the edges connected to this node

	end
end