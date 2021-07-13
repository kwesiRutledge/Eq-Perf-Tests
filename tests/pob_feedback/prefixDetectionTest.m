function tests = prefixDetectionTest
	%Descritption:
	%	Provides tests through the local functions in this file.

	tests = functiontests(localfunctions);

function [sys_out,L1,P_u,P_x0] = get_simple_1dsys()
	%Description:

	L1 = Language([1,2,1,2],[1,1,1,1],[3,1,1,1]);
	T = length(L1.words{1});

	A1 = 0.9;
	B1 = 1;
	C1 = 1; %[1,0];
	
	n_x = size(A1,1);
	n_u = size(B1,2);
	n_y = size(C1,1);

	eta_v = 0.1; eta_w = 0.2;
	Pv1 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));
	Pw1 = Polyhedron('lb',-eta_w*ones(1,n_x) ,'ub',eta_w*ones(1,n_x));
	eta_v = 0.3;
	Pv2 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));

	eta_u = 1.0; eta_x0 = 0.3;
	P_u = Polyhedron('lb',-eta_u*ones(1,n_u) ,'ub',eta_u*ones(1,n_u));
	P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));

	f1 = eta_w;
	f2 = -eta_w;

	aff_dyn_list = [	Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f1,C1,Pw1,Pv2) ];

	sys_out = LCSAS( aff_dyn_list , L1 , 'X0' , P_x0 , 'U' , P_u );

function [BG_out] = get_unconnected_bg1()
	%get_unconnected_bg1
	%Description:
	%	Creates a very particular BeliefGraph to test the BeliefLanguage construction.

	%% Constants %%
	T = 4;
	[lcsas,~,P_u,P_x0] = get_simple_1dsys();

	%% Create Some Fake Consistency Sets
	c_set1 = Polyhedron('lb',[0],'ub',1);
	c_set2 = Polyhedron('lb',-1,'ub',0);

	%% Create Graph %%
	bn1 = BeliefNode(Language([1,2]),0,'ConsistencySet',c_set1);
	bn3 = BeliefNode(Language([1,2]),1);
	bn5 = BeliefNode(Language([1,2]),2);
	bn7 = BeliefNode(Language([1,2]),3);

	bn2 = BeliefNode(Language([2]),0,'ConsistencySet',c_set2);
	bn4 = BeliefNode(Language([2]),1);
	bn6 = BeliefNode(Language([2]),2);
	bn8 = BeliefNode(Language([2]),3);

	edges = [ 	1,3;
				3,5;
				5,7;
				2,4;
				4,6;
				6,8 ];

	BG_out = BeliefGraph(lcsas,P_u,P_x0,'return_empty',true);
	BG_out.N = [bn1,bn2,bn3,bn4,bn5,bn6,bn7,bn8];
	BG_out.E = edges;

function test_prefix_detection1(testCase)
	%Description:

	%% Constants

	BG1 = get_unconnected_bg1();
	controller1 = POB_Feedback(BG1,{1},{2});

	%% Test

	assert( 2 == controller1.prefix_detection(0.1) )