%test_InternalBehaviorSet.m
function tests = internalBehaviorSetsTest
	%disp(localfunctions)
	tests = functiontests(localfunctions);

function include_relevant_libraries()
	%Description:
	%	Attempts to add the relevant libraries/toolboxes to the path.

	%% Constants

	is_on_personal_mac = strcmp(getenv('USER'),'kwesirutledge');
	is_on_great_lakes = (strcmp(getenv('USER'),'krutledg') && isunix) ;

	addpath(genpath('../../'))

	%% Algorithm
	% if is_on_personal_mac
	% 	include_fcns2({'mosek','gurobi','tbxmanager'},'PathToDirectoryWithToolboxes','../../../','Verbose',0)
	% elseif is_on_great_lakes
	% 	include_fcns2({'tbxmanager','YALMIP'},'PathToDirectoryWithToolboxes','../../../','Verbose',0)
	% else
	% 	warning('Warning: There is no guaranteed performance for systems that are not on Kwesi''s computer.')
	% end

function [lcsas,eta_w,eta_v,eta_u] = get_simple_lcsas1()
	%Description:
	%	Retrieves a simple LCSAS in two dimensions for testing.

	T = 4;
	L1 = Language( ones(1,T) , 2*ones(1,T) );

	%Create a simple Language Constrainted Switching System
	A1 = eye(2); A2 = 2*eye(2);
	B1 = [0;1];
	C1 = [1,0];

	n_x = size(A1,1);
	n_u = size(B1,2);
	n_y = size(C1,1);

	eta_v = 0.1; eta_w = 0.2;
	Pv1 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));
	Pw1 = Polyhedron('lb',-eta_w*ones(1,n_x) ,'ub',eta_w*ones(1,n_x));

	eta_u = 0; eta_x0 = 0.3;
	P_u = Polyhedron('lb',-eta_u*ones(1,n_u) ,'ub',eta_u*ones(1,n_u));
	P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));

	f1 = zeros(2,1);

	aff_dyn_list = [	Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1), ...
						Aff_Dyn(A2,B1,f1,C1,Pw1,Pv1) ];

	lcsas = LCSAS( aff_dyn_list , L1 , 'X0' , P_x0 , 'U' , P_u );

function test_InternalBehaviorSet_Constructor1(testCase)
	%test_InternalBehaviorSet_Constructor1.m
	%Description:
	%	This test script is meant to evaluat the constructor for InternalBehaviorSet.

	include_relevant_libraries();

	% Constants
	[lcsas0,~,~,~] = get_simple_lcsas1(); %Get Simple Dynamics
	L1 = Language(lcsas0.L.words{1});
	temp_belief_sequence = [L1;L1];

	%Create Example InternalBehaviorSet
	ibs1 = InternalBehaviorSet(lcsas0,temp_belief_sequence);

	assert( ibs1.t == 2 )

function test_InternalBehaviorSet_Constructor2(testCase)
	%test_InternalBehaviorSet_Constructor2.m
	%Description:
	%	This test script is meant to evaluate the constructor for InternalBehaviorSet()
	%	checking the A matrices.

	include_relevant_libraries();

	% Constants
	[lcsas0,~,~,~] = get_simple_lcsas1(); %Get Simple Dynamics
	L1 = Language(lcsas0.L.words{1});
	temp_belief_sequence = [L1];

	%Create Example InternalBehaviorSet
	ibs1 = InternalBehaviorSet(lcsas0,temp_belief_sequence);

	x0 = 0.2*ones(2,1);
	w0 = 0.1*ones(2,1);
	x1 = lcsas0.Dyn(1).A * x0 + w0;

	u0 = zeros(1,1);

	temp_poly = ibs1.ToPolyhedron();

	assert( temp_poly.contains([x0;x1;u0;w0;x0]) )

function [lcsas,eta_w,eta_v,eta_u] = get_simple_lcsas2()
	%Description:
	%	Retrieves a simple LCSAS in two dimensions for testing.
	%	The two words should lead to an interesting projection set.

	T = 4;
	L1 = Language( ones(1,T) , 2*ones(1,T) );

	%Create a simple Language Constrainted Switching System
	A1 = eye(2);
	B1 = [0;1];
	C1 = [1,0];

	n_x = size(A1,1);
	n_u = size(B1,2);
	n_y = size(C1,1);

	eta_v = 0.1; eta_w = 0.2;
	Pv1 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));
	Pw1 = Polyhedron('lb',-eta_w*ones(1,n_x) ,'ub',eta_w*ones(1,n_x));
	Pw2 = eta_w*[1;0] + Pw1;

	eta_u = 0; eta_x0 = 0.3;
	P_u = Polyhedron('lb',-eta_u*ones(1,n_u) ,'ub',eta_u*ones(1,n_u));
	P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));

	f1 = zeros(2,1);

	aff_dyn_list = [	Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f1,C1,Pw2,Pv1) ];

	lcsas = LCSAS( aff_dyn_list , L1 , 'X0' , P_x0 , 'U' , P_u );

function test_InternalBehaviorSet_ToW1(testCase)
	%test_InternalBehaviorSet_ToW1
	%Description:
	%	This test script is meant to evaluate the method ToW() for InternalBehaviorSet().
	%	In this case, we make sure that it contains the proper hyperrectangle for the simple
	%	lcsas in get_simple_lcsas2().

	include_relevant_libraries();

	% Constants
	[lcsas0,eta_w,~,~] = get_simple_lcsas2(); %Get Simple Dynamics
	L_all = lcsas0.L;
	temp_belief_sequence = [L_all];

	%Create Example InternalBehaviorSet
	ibs1 = InternalBehaviorSet(lcsas0,temp_belief_sequence);

	temp_W = ibs1.ToW();
	W_expected = Polyhedron('lb',[0,-eta_w],'ub',[eta_w,eta_w]);

	assert( (length(temp_W) == 2) && ( W_expected <= temp_W(1) ) && ( W_expected >= temp_W(1) ) )

function test_InternalBehaviorSet_CoversInputPolyhedron1(testCase)
	%Description:
	%	Tests if a given input Polyhedron is covered by a set of InternalBehaviorSets

	%% Constants

	% Get Simple System
	[lcsas0,~,~,~] = get_simple_lcsas1(); %Get Simple Dynamics
	L1 = Language(lcsas0.L.words{1});
	temp_belief_sequence = [L1];

	% Create Target Set to Cover
	W = Polyhedron('lb',[-1,-2],'ub',[1,2]);

	% Create Example IBS
	ibs1 = InternalBehaviorSet(lcsas0,temp_belief_sequence,'ReturnEarly',true,'A',[eye(2);-eye(2)],'b',[1;2;1;0]);
	ibs2 = InternalBehaviorSet(lcsas0,temp_belief_sequence,'ReturnEarly',true,'A',[eye(2);-eye(2)],'b',[1;0;1;2]);

	ibs_array = [ibs1;ibs2];

	assert( ibs_array.CoversInputPolyhedron(W) )

function test_InternalBehaviorSet_CoversInputPolyhedron2(testCase)
	%Description:
	%	Tests if a given input Polyhedron is covered by a set of InternalBehaviorSets

	%% Constants

	% Get Simple System
	[lcsas0,~,~,~] = get_simple_lcsas1(); %Get Simple Dynamics
	L1 = Language(lcsas0.L.words{1});
	temp_belief_sequence = [L1];

	% Create Target Set to Cover
	W = Polyhedron('lb',[-1,-2],'ub',[1,2]);

	% Create Example IBS
	ibs1 = InternalBehaviorSet(lcsas0,temp_belief_sequence,'ReturnEarly',true,'A',[eye(2);-eye(2)],'b',[1;2;1;-1]);
	ibs2 = InternalBehaviorSet(lcsas0,temp_belief_sequence,'ReturnEarly',true,'A',[eye(2);-eye(2)],'b',[1;-1;1;2]);

	ibs_array = [ibs1;ibs2];

	assert( ~ibs_array.CoversInputPolyhedron(W) )

function test_InternalBehaviorSet_CoversInputW1(testCase)
	%Description:
	%	Tests if a given input Polyhedron is covered by a set of InternalBehaviorSets

	%% Constants

	% Get Simple System
	[lcsas0,eta_w,~,~] = get_simple_lcsas2(); %Get Simple Dynamics
	L1 = Language(lcsas0.L.words{1});
	L2 = Language(lcsas0.L.words{2});
	temp_belief_sequence1 = [L1];
	temp_belief_sequence2 = [L2];

	%Create Example InternalBehaviorSet
	ibs1 = InternalBehaviorSet(lcsas0,temp_belief_sequence1);
	ibs2 = InternalBehaviorSet(lcsas0,temp_belief_sequence2);

	% Create Target Set to Cover
	W_prime = Polyhedron('lb',-eta_w*ones(1,2),'ub',eta_w*[2,1]);

	ibs_array = [ibs1;ibs2];

	assert( ibs_array.CoversInputW(W_prime) )

function test_InternalBehaviorSet_CoversInputW2(testCase)
	%Description:
	%	Tests if a given input Polyhedron is covered by a set of InternalBehaviorSets
	%	The input W_prime should NOT BE covered

	%% Constants

	% Get Simple System
	[lcsas0,eta_w,~,~] = get_simple_lcsas2(); %Get Simple Dynamics
	L1 = Language(lcsas0.L.words{1});
	L2 = Language(lcsas0.L.words{2});
	temp_belief_sequence1 = [L1];
	temp_belief_sequence2 = [L2];

	%Create Example InternalBehaviorSet
	ibs1 = InternalBehaviorSet(lcsas0,temp_belief_sequence1);
	ibs2 = InternalBehaviorSet(lcsas0,temp_belief_sequence2);

	% Create Target Set to Cover
	W_prime = Polyhedron('lb',-eta_w*ones(1,2),'ub',eta_w*[2,2]);

	ibs_array = [ibs1;ibs2];

	assert( ~ibs_array.CoversInputW(W_prime) )
