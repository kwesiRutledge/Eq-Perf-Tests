function tests = internal_behavior_set_test
	%disp(localfunctions)
	tests = functiontests(localfunctions);

function test_internal_behavior_set1(testCase)
	%Description:
	%	Tests that the internal behavior set is correctly constructed
	%	when the belief sequence is simple (single word, length one).

	%% Constants

	[lcsas,P_u] = get_1d_lcsas();
	lcsas.U = P_u;


	tempL = Language(lcsas.L.words{1});
	KnowlSequence = repmat(tempL,1,1);

	ApproxIBS1 = lcsas.X0 * Polyhedron('lb',-30,'ub',30) * lcsas.U * lcsas.Dyn(1).P_w * lcsas.X0;

	%% Algorithm 

	[ internal_behavior_set ] = lcsas.internal_behavior_set( KnowlSequence )

	disp([ 'internal_behavior_set.Dim = ' num2str(internal_behavior_set.Dim) ])
	internal_behavior_set.isEmptySet;
	disp([ 'ApproxIBS1.Dim = ' num2str(ApproxIBS1.Dim) ])
	ApproxIBS1.isEmptySet;

	assert( ApproxIBS1 >= internal_behavior_set )

function test_internal_behavior_set2(testCase)
	%Description:
	%	Tests that the internal behavior set is correctly constructed
	%	when the belief sequence is simple (single word, length two).

	%% Constants

	[lcsas,P_u] = get_1d_lcsas();
	lcsas.U = P_u;


	tempL = Language(lcsas.L.words{1});
	KnowlSequence = repmat(tempL,2,1);

	ApproxIBS1 = lcsas.X0 * Polyhedron('lb',-10,'ub',10) * Polyhedron('lb',-30,'ub',30) * lcsas.U * lcsas.U * lcsas.Dyn(1).P_w * lcsas.Dyn(1).P_w * lcsas.X0;

	%% Algorithm 

	[ internal_behavior_set ] = lcsas.internal_behavior_set( KnowlSequence )

	disp([ 'internal_behavior_set.Dim = ' num2str(internal_behavior_set.Dim) ])
	internal_behavior_set.isEmptySet
	disp([ 'ApproxIBS1.Dim = ' num2str(ApproxIBS1.Dim) ])
	ApproxIBS1.isEmptySet

	assert( ApproxIBS1 >= internal_behavior_set )

function test_internal_behavior_set3(testCase)
	%Description:
	%	Tests that the internal behavior set is correctly constructed
	%	when the belief sequence is a bit more compelx (two words, length one).

	%% Constants

	[lcsas,P_u] = get_1d_lcsas();
	lcsas.U = P_u;


	tempL = Language(lcsas.L.words{[1:2]});
	KnowlSequence = repmat(tempL,1,1);

	ApproxIBS1 = lcsas.X0 * Polyhedron('lb',-10,'ub',10) * lcsas.U * lcsas.Dyn(1).P_w * lcsas.Dyn(1).P_w * lcsas.X0 * lcsas.X0;

	%% Algorithm 

	[ internal_behavior_set ] = lcsas.internal_behavior_set( KnowlSequence );

	disp([ 'internal_behavior_set.Dim = ' num2str(internal_behavior_set.Dim) ])
	internal_behavior_set.isEmptySet
	disp([ 'ApproxIBS1.Dim = ' num2str(ApproxIBS1.Dim) ])
	ApproxIBS1.isEmptySet

	assert( ApproxIBS1 >= internal_behavior_set )

function test_internal_behavior_set4(testCase)
	%Description:
	%	Tests that the internal behavior set is correctly constructed
	%	when the belief sequence is complex and should lead to empty set.

	%% Constants

	[lcsas,P_u] = get_1d_lcsas();
	lcsas.Dyn(2).P_w = Polyhedron('lb',-0.5,'ub',-0.1);
	lcsas.U = P_u;


	tempL = Language(lcsas.L.words{[1:2]});
	KnowlSequence = repmat(tempL,2,1);

	%% Algorithm 

	[ internal_behavior_set ] = lcsas.internal_behavior_set( KnowlSequence );

	disp([ 'internal_behavior_set.Dim = ' num2str(internal_behavior_set.Dim) ])
	internal_behavior_set.isEmptySet

	assert( internal_behavior_set.isEmptySet )