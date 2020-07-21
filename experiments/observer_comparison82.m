function [results] = observer_comparison81( varargin )
	%observer_comparison82.m
	%Description:
	%	Playing around with the pendulum + elastic wall example from Tobia Marcucci's paper.
	%

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	test_name = 'observer_comparison82';

	disp(['Beginning ' test_name '.' ])
	disp('Testing Pendulum with Elastic Walls')
	disp(' ')

	disp('Defining Constants.')

	m = 1;
	l = 1;
	g = 10;
	k = 100;
	d = .1;
	h = .01;

	%% Create Dynamics Matrices
	%
	%	x-dot = { A1 x + B1 u, 		if (x,u) \in D1
	%			{ A2 x + B2 u + c2,	if (x,u) \in D2
	%

	A1 = [ 	0 , 1 ;
			g/l , 0 ];

	B1 = [ 0 ; 1/(m*l^2) ];

	A2 = [ 	0 , 1 ;
			(g/l) - (k/m), 0 ];

	B2 = B1;

	c2 = [0; (k*d)/(m*l)];

	%% Create the Domains of each of these functions
	n = 2;
	eta_x = 1.5;
	X_template = Polyhedron('lb',-eta_x*ones(1,n),'ub',eta_x*ones(1,n));
	eta_u = 4;
	Pu = Polyhedron('lb',-eta_u*ones(1,n),'ub',eta_u*ones(1,n));

	D1 = X_template.intersect( Polyhedron('A',[1,0],'b',d/l) );
	D1 = D1 * Pu;

	D2 = X_template.intersect( Polyhedron('A',-[1,0],'b',-d/l) );
	D2 = D2 * Pu;


	eta_w = 0.4;
	Pw_template = Polyhedron('lb',-eta_w*ones(1,n),'ub',eta_w*ones(1,n));
	Pw1 = Pw_template + [eta_w;0];
	Pw2 = Pw_template + [0;eta_w];
	Pv = Polyhedron('lb',zeros(1,n),'ub',zeros(1,n));

	ad1 = Aff_Dyn(A,B,f,C,Pw1,Pv);
	ad2 = Aff_Dyn(A,B,f,C,Pw2,Pv);

	sys = LCSAS( [ad1,ad2] , Language([1,1,1],[2,2,2]) );

	eta_u = 1;
	Pu = Polyhedron('lb',-eta_u*ones(1,m),'ub',eta_u*ones(1,m));

	eta_x0 = 0.1;
	Px0 = Polyhedron('lb',-eta_x0*ones(1,n),'ub',eta_x0*ones(1,n));

	results.sys = sys;

	verbosity = 1;

	flags.ConstructBeliefGraph = false;

	%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Display Belief Graph %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%

	if flags.ConstructBeliefGraph
		%% Get the Belief Graph for this Simple 2D System
		BG = BeliefGraph(sys,Pu,Px0,'verbosity',verbosity);

		figure;
		BG.plot()
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Optimize Using Binary Variables to Incorporate Reachability/Consistency Results %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	disp('Beginning test 1.')
	disp('Creating an optimization that uses the mixed integer linear variables.')
	

end