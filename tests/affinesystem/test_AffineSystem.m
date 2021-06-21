%test_AffineSystem.m
%Description:
%	Tests the constructor for AffineSystem.

function tests = test_AffineSystem
	%disp(localfunctions)
	tests = functiontests(localfunctions);
end

function include_relevant_libraries()
	%Description:
	%	Attempts to add the relevant libraries/toolboxes to the path.

	%% Constants

	%% Algorithm
	addpath(genpath('../../functions'));
end

function test1_Constructor(testCase)

	% Include Relevant Libraries
	include_relevant_libraries()

	A = 1;
	B = 3;
	f = 2;

	as1 = AffineSystem(A,B,f);

	assert( (as1.A == A) && (as1.B == B ) && (as1.f == f) && ...
			(as1.C == 0) )

end