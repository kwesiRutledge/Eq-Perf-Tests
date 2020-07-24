function tests = testBeliefNode
	%disp(localfunctions)
	tests = functiontests(localfunctions);

function test1_BeliefNode(testCase)
	%Description:
	%	Testing the empty constructor.

	L1 = Language([1],[2],[3]);

	bn1 = BeliefNode( L1 , 1 );

	assert( (L1 == bn1.subL) & (bn1.t == 1) & isempty(bn1.c_set) & isempty(bn1.FullTrajectorySet) )

function test2_BeliefNode(testCase)
	%Description:
	%	Testing the novel constructors.

	L1 = Language([1],[2],[3]);

	eta1 = 0.1;
	P1 = Polyhedron('lb',-eta1*ones(1,2),'ub',eta1*ones(1,2));

	bn1 = BeliefNode( L1 , 1 , P1 );

	bn1.c_set;

	assert( (L1 == bn1.subL) & (bn1.t == 1) & (P1 <= (bn1.c_set) ) & (P1 <= (bn1.c_set)) )

function test3_BeliefNode_eq(testCase)
	%Description:
	%	Tests the "==" operator for Belief Nodes.

	L1 = Language([1],[2],[3]);
	L2 = Language([2],[3]);	
	L3 = Language([1],[2],[3]);

	bn1 = BeliefNode( L1 , 0 );
	bn2 = BeliefNode( L2 , 0 );
	bn3 = BeliefNode( L3 , 0 );
	bn4 = BeliefNode( L1 , 1 );

	assert( (bn1 == bn3) & ~(bn1 == bn2) & ~(bn1 == bn4) )
	