function tests = testBeliefNode
	disp(localfunctions)
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

	assert( (L1 == bn1.subL) & (bn1.t == 1) & P1 == bn1.c_set )

	