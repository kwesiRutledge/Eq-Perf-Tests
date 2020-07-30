function tests = testBeliefNode
	%disp(localfunctions)
	tests = functiontests(localfunctions);

function test1_BeliefNode(testCase)
	%Description:
	%	Testing the empty constructor.

	%% Including Libraries
	include_relevant_libraries();

	%% Algorithm

	L1 = Language([1],[2],[3]);

	bn1 = BeliefNode( L1 , 1 );

	assert( (L1 == bn1.subL) & (bn1.t == 1) & isempty(bn1.c_set) & isempty(bn1.FullTrajectorySet) )

function test2_BeliefNode(testCase)
	%Description:
	%	Testing the novel constructors.

	%% Including Libraries
	include_relevant_libraries();

	%% Algorithm

	L1 = Language([1],[2],[3]);

	eta1 = 0.1;
	P1 = Polyhedron('lb',-eta1*ones(1,2),'ub',eta1*ones(1,2));

	bn1 = BeliefNode( L1 , 1 , P1 );

	bn1.c_set;

	assert( (L1 == bn1.subL) & (bn1.t == 1) & (P1 <= (bn1.c_set) ) & (P1 <= (bn1.c_set)) )

function test3_BeliefNode_eq(testCase)
	%Description:
	%	Tests the "==" operator for Belief Nodes.

	%% Including Libraries
	include_relevant_libraries();

	%% Algorithm

	L1 = Language([1],[2],[3]);
	L2 = Language([2],[3]);	
	L3 = Language([1],[2],[3]);

	bn1 = BeliefNode( L1 , 0 );
	bn2 = BeliefNode( L2 , 0 );
	bn3 = BeliefNode( L3 , 0 );
	bn4 = BeliefNode( L1 , 1 );

	assert( (bn1 == bn3) & ~(bn1 == bn2) & ~(bn1 == bn4) )
	
function include_relevant_libraries()
	%Description:
	%	Attempts to add the relevant libraries/toolboxes to the path.

	%% Constants

	is_on_personal_mac = strcmp(getenv('USER'),'kwesirutledge');
	is_on_great_lakes = (strcmp(getenv('USER'),'krutledg') && isunix) ;

	addpath('../')

	verbose = 0;

	%% Algorithm
	if is_on_personal_mac
		include_fcns2({'mosek','gurobi','tbxmanager'},'PathToDirectoryWithToolboxes','../../../','Verbose',verbose)
	elseif is_on_great_lakes
		include_fcns2({'tbxmanager','YALMIP'},'PathToDirectoryWithToolboxes','../../../')
	else
		warning('Warning: There is no guaranteed performance for systems that are not on Kwesi''s computer.')
	end

	addpath(genpath('../functions'));

function test1_get_consistency_causing_disturbances(testCase)
	%Description:
	%	This test verifies that the error handling of get_consistency_causing_disturbances().
	%

	%% Including Libraries
	include_relevant_libraries();

	%% Constants
	L1 = Language([1],[2],[3]);
	bn1 = BeliefNode( L1 , 0 );

	%% Algorithms
	try
		bn1.get_consistency_causing_disturbances()
		assert(false)
	catch e
		assert(strcmp(e.message,'This function cannot be called when the FullTrajectorySet is not defined for this BeliefNode.'));
	end

function test2_get_consistency_causing_disturbances(testCase)
	%Description:
	%	This test verifies get_consistency_causing_disturbances() retrieves
	%	the proper polyhedron when the node's subL.cardinality() = 1.
	%

	%% Including Libraries
	include_relevant_libraries();

	%% Constants
	[lcsas2,P_u,Pw1,Pw2,eta_v,eta_x0] = get_1d_lcsas();

	BG = BeliefGraph(lcsas2,P_u,lcsas2.X0,'return_empty',true);
	initial_node = BG.get_initial_beliefnodes();
	% disp('initial_node.subL = ')
	% disp(initial_node.subL.words)

	level1 = BG.post( initial_node , P_u , lcsas2.X0 );

	% disp(['length(level1) = ' num2str(length(level1))])
	% disp(['level1(1).subL.cardinality() = ' num2str(level1(1).subL.cardinality()) ])
	% disp(['level1(1).subL.words{1} = ' num2str(level1(1).subL.words{1}) ])
	% disp(['level1(1).c_set.Dim = ' num2str(level1(1).c_set.Dim) ])
	

	bn2 = level1(1);

	% superset_of_c_set = Polyhedron('A', [ 	1  , 0  , 0 ;
	% 									-1 , 0  , 0 ;
	% 									0  , 1  , 0 ;
	% 									0  , -1 , 0 ], ...
	% 							'b', [ eta_x0 ; eta_x0 ; P_u.b ;  ])

	superset_of_c_set = Polyhedron('lb',[-eta_x0-eta_v,0.9*(-eta_x0)-P_u.b(1)-eta_v-Pw2.b(1),-P_u.b(1)], ...
	 							'ub',[eta_x0+eta_v,0.9*(eta_x0)+P_u.b(2)+eta_v+Pw2.b(2),P_u.b(1)])

	expected_causing_disturbances = ...
		Pw2 * Polyhedron('lb',-eta_v,'ub',eta_v) * lcsas2.X0;

	% disp(bn2.c_set.b)
	% disp(bn2.c_set.A)

	%% Algorithms
	result2 = bn2.get_consistency_causing_disturbances(lcsas2);

	% disp( superset_of_c_set <= bn2.c_set )
	% disp( superset_of_c_set >= bn2.c_set )

	% disp(result2 <= expected_causing_disturbances)
	% disp(result2 >= expected_causing_disturbances)

	assert( (bn2.c_set <= superset_of_c_set) && ...
			(result2 == expected_causing_disturbances) )

function [lcsas_out,P_u,Pw1,Pw2,eta_v,eta_x0] = get_1d_lcsas_v2( varargin )
	%Description:
	%	Create a one-dimensional system to quickly synthesize belief trees.
	%
	%Usage:
	%	[lcsas_out,P_u] = get_1d_lcsas()

	%% Constants %%

	a = 0.9;
	b = 1;

	Pw1 = Polyhedron('lb',0.0,'ub',0.5);
	Pw2 = Polyhedron('lb',-0.5,'ub',0.0);

	eta_v = 0.25;
	Pv = Polyhedron('lb',-eta_v,'ub',eta_v);

	eta_x0 = 0.25;
	P_x0 = Polyhedron('lb',-eta_x0,'ub',eta_x0);

	%% Algorithm %%

	eta_u = 1;
	P_u = Polyhedron('lb',-eta_u,'ub',eta_u);

	ad1 = Aff_Dyn(a,b,zeros(1),eye(1),Pw1,Pv);
	ad2 = Aff_Dyn(a,b,zeros(1),eye(1),Pw2,Pv);

	lcsas_out = LCSAS( [ad1,ad2] , Language([1,1,1],[2,2,2],[1,2,1]) , 'X0' , P_x0 );

function test3_get_consistency_causing_disturbances(testCase)
	%Description:
	%	This test verifies that the error handling of get_consistency_causing_disturbances().
	%

	%% Including Libraries
	include_relevant_libraries();

	%% Constants
	[lcsas3,P_u,Pw1,Pw2,eta_v,eta_x0] = get_1d_lcsas_v2();

	BG = BeliefGraph(lcsas3,P_u,lcsas3.X0,'return_empty',true);
	initial_node = BG.get_initial_beliefnodes();
	disp('initial_node.subL = ')
	disp(initial_node.subL.words)

	level1 = BG.post( initial_node , P_u , lcsas3.X0 );

	disp(['length(level1) = ' num2str(length(level1))])
	disp(['level1(2).subL.cardinality() = ' num2str(level1(2).subL.cardinality()) ])
	disp(['level1(2).subL.words{1} = ' num2str(level1(2).subL.words{1}) ])
	disp(['level1(2).subL.words{1} = ' num2str(level1(2).subL.words{2}) ])
	disp(['level1(2).c_set.Dim = ' num2str(level1(2).c_set.Dim) ])
	

	bn3 = level1(2);

	expected_causing_disturbances = ...
		Pw1 * Polyhedron('lb',-eta_v,'ub',eta_v) * lcsas3.X0;

	disp(bn3.c_set.b)
	disp(bn3.c_set.A)

	%% Algorithms
	result3 = bn3.get_consistency_causing_disturbances(lcsas3)

	for result_idx = 1:length(result3)
		disp(result3(result_idx) <= expected_causing_disturbances)
		disp(result3(result_idx) >= expected_causing_disturbances)
	end


	assert( ( length(result3) == bn3.subL.cardinality() ) && ...
			(result3(1) <= expected_causing_disturbances) )