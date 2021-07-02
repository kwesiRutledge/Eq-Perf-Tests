function tests = testConsistencySet
	%disp(localfunctions)
	tests = functiontests(localfunctions);

function test1_ConsistencySet(testCase)
	%Description:
	%	

	L1 = Language([1,2,1,2],[3,4,3,4],[5,1,1,1]);
	T = length(L1.words{1});

	A1 = [0,1;0.1,-0.05];
	B1 = [0;1];
	C1 = [1,0];
	
	n_x = size(A1,1);
	n_u = size(B1,2);
	n_y = size(C1,1);

	eta_v = 0.1; eta_w = 0.2;
	Pv1 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));
	Pw1 = Polyhedron('lb',-eta_w*ones(1,n_x) ,'ub',eta_w*ones(1,n_x));
	eta_v = 0.3;
	Pv2 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));

	eta_u = 0; eta_x0 = 0.3;
	P_u = Polyhedron('lb',-eta_u*ones(1,n_u) ,'ub',eta_u*ones(1,n_u));
	P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));

	f1 = eta_w*[zeros(n_x-1,1);1];
	f2 = eta_w*[1;zeros(n_x-1,1)];
	f3 = -f1;
	f4 = -f2;

	aff_dyn_list = [	Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f3,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f4,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f1,C1,Pw1,Pv2) ];

	lcsas0 = LCSAS( aff_dyn_list , L1 , 'X0' , P_x0 );

	try
		[ConsistencySet1] = ExternalBehaviorSet( lcsas0 , L1 , 'fb_method' , 'output' );
	catch e
		assert(strcmp(e.message,'The field U is empty. Please define it as a Polyhedron.'))
	end

function test3_ConsistencySet(testCase)
	%Description:
	%	

	L1 = Language([1,2,1,2],[3,4,3,4],[5,1,1,1]);
	T = length(L1.words{1});

	A1 = [0,1;0.1,-0.05];
	B1 = [0;1];
	C1 = [1,0];
	
	n_x = size(A1,1);
	n_u = size(B1,2);
	n_y = size(C1,1);

	eta_v = 0.1; eta_w = 0.2;
	Pv1 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));
	Pw1 = Polyhedron('lb',-eta_w*ones(1,n_x) ,'ub',eta_w*ones(1,n_x));
	eta_v = 0.3;
	Pv2 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));

	eta_u = 0; eta_x0 = 0.3;
	P_u = Polyhedron('lb',-eta_u*ones(1,n_u) ,'ub',eta_u*ones(1,n_u));
	P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));

	f1 = eta_w*[zeros(n_x-1,1);1];
	f2 = eta_w*[1;zeros(n_x-1,1)];
	f3 = -f1;
	f4 = -f2;

	aff_dyn_list = [	Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f3,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f4,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f1,C1,Pw1,Pv2) ];

	lcsas0 = LCSAS( aff_dyn_list , L1 , 'X0' , P_x0 , 'U' , P_u );

	[ConsistencySet1] = ExternalBehaviorSet( lcsas0 , L1 , 'fb_method' , 'output' );
	[ConsistencySet2] = ExternalBehaviorSet( lcsas0 , L1 , 'fb_method' , 'output' );

	assert( (ConsistencySet1 <= ConsistencySet2) & (ConsistencySet1 >= ConsistencySet2) )

function test2_ConsistencySet(testCase)
	%Description:
	%	Testing the conversion of consistency sets to times.

	L1 = Language([1,2,1,2],[3,4,3,4],[5,1,1,1]);
	T = length(L1.words{1});

	A1 = [0,1;0.1,-0.05];
	B1 = [0;1];
	C1 = [1,0];
	
	n_x = size(A1,1);
	n_u = size(B1,2);
	n_y = size(C1,1);

	eta_v = 0.1; eta_w = 0.2;
	Pv1 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));
	Pw1 = Polyhedron('lb',-eta_w*ones(1,n_x) ,'ub',eta_w*ones(1,n_x));
	eta_v = 0.3;
	Pv2 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));

	eta_u = 0; eta_x0 = 0.3;
	P_u = Polyhedron('lb',-eta_u*ones(1,n_u) ,'ub',eta_u*ones(1,n_u));
	P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));

	f1 = eta_w*[zeros(n_x-1,1);1];
	f2 = eta_w*[1;zeros(n_x-1,1)];
	f3 = -f1;
	f4 = -f2;

	aff_dyn_list = [	Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f3,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f4,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f1,C1,Pw1,Pv2) ];

	lcsas0 = LCSAS( aff_dyn_list , L1 , 'X0' , P_x0 );

	temp_L = Language(L1.words{1});

	[ConsistencySet1,InternalBehaviorSet1] = lcsas0.consistent_set( 1 , temp_L , P_u , P_x0 , 'fb_method' , 'output' );
	[ConsistencySet2,InternalBehaviorSet2] = lcsas0.consistency_set2( 1 , temp_L , P_u , P_x0 , 'fb_method' , 'output' );

	found_t_fromv1_cset = (1 == lcsas0.behavior_set2time( ConsistencySet1 , 'ConsistencySet' , 'ConsistencySetVersion' , 1 , 'FeedbackMethod' , 'output' ));
	found_t_fromv1_ebset = (1 == lcsas0.behavior_set2time( InternalBehaviorSet1 , 'InternalBehaviorSet_1word' , 'ConsistencySetVersion' , 1 , 'FeedbackMethod' , 'output' ));
	found_t_fromv2_cset = (1 == lcsas0.behavior_set2time( ConsistencySet2 , 'ConsistencySet' , 'ConsistencySetVersion' , 2 , 'FeedbackMethod' , 'output' ));
	found_t_fromv2_ebset = (1 == lcsas0.behavior_set2time( InternalBehaviorSet2 , 'InternalBehaviorSet_1word' , 'ConsistencySetVersion' , 2 , 'FeedbackMethod' , 'output' ));

	assert( found_t_fromv1_cset & found_t_fromv1_ebset & found_t_fromv2_cset & found_t_fromv2_ebset )
	