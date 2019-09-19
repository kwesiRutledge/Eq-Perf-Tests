function [results] = observer_comparison58( varargin )
	%Description:
	%	Testing the all mode observer construction algorithm

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	dim = 2;

	L1 = Language([1,1,2,3],[1,2,2,4],[1,2,4,3]);

	%Create a simple Language Constrainted Switching System
	A1 = eye(dim);
	B1 = eye(dim);
	C1 = eye(dim);
	f1 = [0;1];

	eta_v = 0; eta_w = 0.5;
	Pv1 = Polyhedron('lb',-eta_v*ones(1,dim) ,'ub',eta_v*ones(1,dim));
	Pw1 = Polyhedron('lb',-eta_w*ones(1,dim) ,'ub',eta_w*ones(1,dim));

	eta_u = 0.1;
	P_u = Polyhedron('lb',-eta_u*ones(1,dim) ,'ub',eta_u*ones(1,dim));

	ad1 = Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1);

	f2 = [1;0];
	f3 = -f1;
	f4 = -f2;

	lcss = [	ad1,...
				Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1),...
				Aff_Dyn(A1,B1,f3,C1,Pw1,Pv1),...
				Aff_Dyn(A1,B1,f4,C1,Pw1,Pv1)];

	%Define Sets
	eta_u = 0; eta_x0 = 0.3;
	P_u = Polyhedron('lb',-eta_u*ones(1,dim) ,'ub',eta_u*ones(1,dim));
	P_x0 = Polyhedron('lb',-eta_x0*ones(1,dim),'ub',eta_x0*ones(1,dim));

	% Select Matrix
	RT = @(n,t,T) [zeros(n,n*t) eye(n) zeros(n,n*(T-t))];

	%% Synthesizing a feedback gain that tries to recover %%

	% disp('Experiment 1: Creating LCSAS')
	lcsas1 = LCSAS(lcss,L1);
	disp('- Created LCSAS object.')
	disp(' ')

	temp_rot = rotz(45);

	A1 = temp_rot(1:2,1:2);
	B1 = diag([1,2]);
	C1 = eye(dim);

	f = zeros(dim,1);

	temp_rot = rotz(-45);
	A2 = temp_rot(1:2,1:2);
	B2 = diag([2,1]);
	C2 = eye(dim);

	eta_v = 0; eta_w = 0.7;
	Pv2 = Polyhedron('lb',-eta_v*ones(1,dim) ,'ub',eta_v*ones(1,dim));
	Pw2 = Polyhedron('lb',-eta_w*ones(1,dim) ,'ub',eta_w*ones(1,dim));

	ad1 = Aff_Dyn(A1,B1,f,C1,Pw1,Pv1);
	ad2 = Aff_Dyn(A2,B2,f,C2,Pw2,Pv2);

	lcsas2 = LCSAS([ad1,ad2],Language([1,1,2,2],[1,2,1,2]));

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Experiment 1: Testing the All Mode Matrices %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	L1 = place(ad1.A',ad1.C',[0.25,0.1])';
	L2 = place(ad2.A',ad2.C',[0.25,0.1])';

	obsv_gains = L1;
	obsv_gains(:,:,2) = L2;

	disp('Experiment 1: Testing the All Mode Matrices')
	
	[ A_Lt , B_Lt , k_Lt , L_Lt , C_Lt ] = lcsas2.get_allmode_obsv_update_matrices(1,obsv_gains)
	
	disp(' - Created Matrices')
	disp(' ')

	results.exp1.sys = lcsas2;
	results.exp1.A_Lt = A_Lt;
	results.exp1.B_Lt = B_Lt;
	results.exp1.k_Lt = k_Lt;
	results.exp1.L_Lt = L_Lt;
	results.exp1.C_Lt = C_Lt;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Experiment 2: Testing the All Mode MPC Matrices %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	disp('Experiment 1: Testing the All Mode Matrices')
	
	[ S_obsv , H_obsv , Cbar_obsv , J_obsv , k_obsv ] = lcsas2.get_allmode_obsv_mpc_matrices(obsv_gains)
	
	disp(' - Created Matrices')
	disp(' ')

	results.exp1.sys = lcsas2;
	results.exp1.A_Lt = A_Lt;
	results.exp1.B_Lt = B_Lt;
	results.exp1.k_Lt = k_Lt;
	results.exp1.L_Lt = L_Lt;
	results.exp1.C_Lt = C_Lt;

	
end