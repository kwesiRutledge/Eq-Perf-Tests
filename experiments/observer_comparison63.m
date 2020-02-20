function [results] = observer_comparison63( varargin )
	%observer_comparison63.m
	%Description:
	%	Testing the generation of a cover using the modified post algorithm.

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%



	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	dt = 0.1;

	%% Create Continuous Time System Matrices %%

	n = 4;
	m = 2;

	A = zeros(n);
	A(1,2) = 1;
	A(3,4) = 1;

	B = [eye(m);zeros(m)];

	B_w = [zeros(m);eye(m)];

	C = eye(n);
	C_v = eye(n);

	local_eta = 0.1;
	tracking_eta = 0.2;
	P_v = Polyhedron('lb', -[local_eta*ones(1,2),tracking_eta*ones(1,2)] , 'ub' , [local_eta*ones(1,2),tracking_eta*ones(1,2)] );

	human_eta = 0.2;
	P_w = Polyhedron('lb',-human_eta*ones(1,2));

	P_w1 = [0;0.2] + P_w;
	P_w2 = [0;-0.2] + P_w;
	P_w3 = [0.2;0] + P_w;
	P_w4 = [-0.2;0] + P_w;

	%% Create Discretized Versions of the A,B, and B_w matrices %%
	sys1 = ss(A,B,eye(n),0);
	dsys1 = c2d(sys1,dt);

	sys2 = ss(A,B_w,eye(n),0);
	dsys2 = c2d(sys2,dt);

	A_d = dsys1.A;
	B_d = dsys1.B;
	Bw_d = dsys2.B;

	%% Create the set of Aff_Dyn objects %%

	ad1 = Aff_Dyn(	A_d,B_d,zeros(n,1),C,...
					P_w1,P_v, ...
					Bw_d,eye(n));

	ad2 = Aff_Dyn(  A_d,B_d,zeros(n,1),C,...
					P_w2,P_v, ...
					Bw_d,eye(n));

	ad3 = Aff_Dyn(  A_d,B_d,zeros(n,1),C,...
					P_w3,P_v, ...
					Bw_d,eye(n));

	ad4 = Aff_Dyn(  A_d,B_d,zeros(n,1),C,...
					P_w4,P_v, ...
					Bw_d,eye(n));

	%% Create Language for the Switching System %%
	T = 10;
	L = Language([ones(1,T)],4*ones(1,T))

	human_passing_sys = LCSAS([ad1,ad2,ad3,ad4],L);

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	safe_separation = 0.2;

	%% Create the Safe Set that represents when the robot is at least 'safe_separation' away from the human %%
	safe_set1 = Polyhedron('A',-[1,0,-1,0],'b',-safe_separation);
	safe_set2 = Polyhedron('A',[1,0,-1,0],'b',-safe_separation);
	safe_set3 = Polyhedron('A',-[0,1,0,-1],'b',-safe_separation);
	safe_set4 = Polyhedron('A',[0,1,0,-1],'b',-safe_separation);

	safe_set_full = PolyUnion( [safe_set1 , safe_set2 , safe_set3 , safe_set4] );

	%%%%%%%%%%%%%
	%% Results %%
	%%%%%%%%%%%%%
	
	results.System = human_passing_sys;
	results.SafeSet = safe_set_full;
end