function tests = test_p_subseteq_pu
	%disp(localfunctions)
	tests = functiontests(localfunctions);

function test_p_subseteq_pu1(testCase)
	%Description:
	%	Tests that the function pa_subseteq() correctly determines that
	%	a given Polyhedron IS NOT a subset of an array of Polyhedron

	% Constants
	a = Polyhedron('lb',0,'ub',1);
	b = Polyhedron('lb',2,'ub',3);

	c = Polyhedron('lb',1,'ub',2);

	%% Algorithm

	tf = p_subseteq_pu( c , PolyUnion([a,b]) );

	assert(tf == false)

function test_p_subseteq_pu2(testCase)
	%Description:
	%	Tests that the function pa_subseteq() correctly determines that
	%	a given Polyhedron IS a subset of an array of Polyhedron

	% Constants
	a = Polyhedron('lb',0,'ub',2);
	b = Polyhedron('lb',1,'ub',3);

	c = Polyhedron('lb',1,'ub',2);

	%% Algorithm

	tf = p_subseteq_pu( c , PolyUnion([a,b]) );

	assert(tf == true)

function test_p_subseteq_pu3(testCase)
	%Description:
	%	Tests that the function pa_subseteq() correctly determines that
	%	a given Polyhedron IS a subset of an array of Polyhedron

	% Constants
	a = Polyhedron('lb',0,'ub',2);
	b = Polyhedron('lb',1,'ub',3);
	d = Polyhedron('lb',5,'ub',100);

	c = Polyhedron('lb',1,'ub',2);

	%% Algorithm

	tf = p_subseteq_pu( c , PolyUnion([a,b,d]) );

	assert(tf == true)


function test_p_subseteq_pu4(testCase)
	%Description:
	%	Tests that the function pa_subseteq() correctly determines that
	%	a given Polyhedron IS a subset of an array of Polyhedron

	% Constants
	a = Polyhedron('lb',0,'ub',2);
	b = Polyhedron('lb',1,'ub',3);
	d = Polyhedron('lb',5,'ub',100);

	c = Polyhedron('lb',0.75,'ub',2.25);

	%% Algorithm

	tf = p_subseteq_pu( c , PolyUnion([a,b,d]) );

	assert(tf == true)