function [results] = observer_comparison83( varargin )
	%observer_comparison83.m
	%Description:
	%	Constructing Belief Graph for illustration of one of the points
	%	about control synthesis as a polynomial optimization problem.
	%

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	test_name = 'observer_comparison83';

	disp(['Beginning ' test_name '.' ])
	disp('Constructing a belief graph to use for testing polynomial optimization ideas.')
	disp(' ')

	disp('- Defining Constants.')

    %Create LCSLS
    %
    %	x-dot = { A1 x + B1 u + d, 	\sigma=1
	%			{ A2 x + B2 u + d,	\sigma=2
    
    eta_u = 0.5;
    eta_d = 0.25;
    
    %System 1
	A1 = 1;
    B1 = 1;
    D1 = Polyhedron('lb',-eta_d,'ub',eta_d);
    U1 = Polyhedron('lb',-eta_u,'ub',eta_u);
    
    %System 2
    A2 = 0.5;
    B2 = 1;
    D2 = D1;
    U2 = U1;
    
    n = size(A1,1);
    Pv = Polyhedron('lb',0,'ub',0); %Want for this to be state feedback

	ad1 = Aff_Dyn(A1,B1,zeros(n,1),eye(n),D1,Pv);
	ad2 = Aff_Dyn(A2,B2,zeros(n,1),eye(n),D2,Pv);

    T = 4;
	sys = LCSAS( [ad1,ad2] , Language(1*ones(1,T),2*ones(1,T) ) );

	results.Parameters.sys = sys;
    results.Parameters.U1 = U1;
    results.Parameters.U2 = U2;
    
    %Define Initial State Set
    Px0 = Polyhedron('lb',-0.3,'ub',0.3);
    results.Parameters.Px0 = Px0;
    
	%Define Target
	P_target = Polyhedron('lb',1,'ub',2);
	results.Parameters.P_target = P_target;

	% Plotting Information
	verbosity = 1;

	flags.ConstructBeliefGraph = true;

	results.Parameters.verbosity = verbosity;
	results.Parameters.flags = flags;

	%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Display Belief Graph %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%

	if flags.ConstructBeliefGraph
		%% Get the Belief Graph for this Simple 2D System
		BG = BeliefGraph(sys,U1,Px0,'verbosity',verbosity,'fb_method','state');
		results.Experiment0.BeliefGraph = BG;

		disp('- Created BeliefGraph object.')

		figure;
		BG.plot()
	end

	%%%%%%%%%%%%%%%%%
	%% T=1 Problen %%
	%%%%%%%%%%%%%%%%%

	

end