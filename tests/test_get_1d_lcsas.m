%test_get_1d_lcsas
%Description:
%	Tests the simple function get_1d_lcsas()

function tests = test_get_1d_lcsas
	%disp(localfunctions)
	tests = functiontests(localfunctions);
end

function include_relevant_libraries()
	%Description:
	%	Attempts to add the relevant libraries/toolboxes to the path.

	%% Constants

	is_on_personal_mac = strcmp(getenv('USER'),'kwesirutledge');
	is_on_great_lakes = (strcmp(getenv('USER'),'krutledg') && isunix) ;

	addpath('../')

	%% Algorithm
	if is_on_personal_mac
		include_fcns('mosek','gurobi','tbxmanager')
	elseif is_on_great_lakes
		include_fcns('tbxmanager','YALMIP')
	else
		warning('Warning: There is no guaranteed performance for systems that are not on Kwesi''s computer.')
	end

	addpath(genpath('../functions'));
end


function test1_constructor(testCase)
	%Description:
	%

	%% Including Libraries
	include_relevant_libraries();

	%% Constants
	sys1 = get_1d_lcsas();

	Pw1_expected = Polyhedron('lb',0.0,'ub',0.5);
	Pw2_expected = Polyhedron('lb',-0.5,'ub',0.0);

	Px0_expected = Polyhedron('lb',-0.25,'ub',0.25);

	%% Algorithm
	assert(	(sys1.Dyn(1).A == 0.9) && (sys1.Dyn(1).B == 1) && ...
			(sys1.Dyn(1).P_w == Pw1_expected) && ...
			(sys1.Dyn(2).P_w == Pw2_expected) && ...
			(sys1.X0 == Px0_expected) )

end