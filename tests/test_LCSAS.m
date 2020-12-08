%test_LCSAS.m
%Description:
%	Tests the function LCSAS.

function tests = test_LCSAS
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
		include_fcns2({'mosek','gurobi','tbxmanager'},'PathToDirectoryWithToolboxes','../../../','Verbose',0)
	elseif is_on_great_lakes
		include_fcns2({'tbxmanager','YALMIP'},'PathToDirectoryWithToolboxes','../../../','Verbose',0)
	else
		warning('Warning: There is no guaranteed performance for systems that are not on Kwesi''s computer.')
	end

end

function test1_Dim_x(testCase)
	%Description:
	%	Verifies if Dim_x correctly identifies that a one-dimensional sytem is one-dimensional.

	%% Including Libraries
	include_relevant_libraries();

	%% Constants
	sys1 = get_1d_lcsas();

	%% Algorithm
	assert(	sys1.Dim_x() == 1 )

end

function test2_Dim_x(testCase)
	%Description:
	%	Verifies if Dim_x correctly identifies that the spacecraft system is six-dimensional.

	%% Including Libraries
	include_relevant_libraries();

	%% Constants
	sys2 = get_spacecraft_lcsas1();

	%% Algorithm
	assert(	sys2.Dim_x() == 6 )

end

function test1_Dim_u(testCase)
	%Description:
	%	Verifies if Dim_u correctly identifies that a one-dimensional sytem's input is one-dimensional.

	%% Including Libraries
	include_relevant_libraries();

	%% Constants
	sys1 = get_1d_lcsas();

	%% Algorithm
	assert(	sys1.Dim_u() == 1 )

end

function test2_Dim_u(testCase)
	%Description:
	%	Verifies if Dim_u correctly identifies that the spacecraft system's input is two-dimensional.

	%% Including Libraries
	include_relevant_libraries();

	%% Constants
	sys2 = get_spacecraft_lcsas1();

	%% Algorithm
	assert(	sys2.Dim_u() == 2 )

end

function test1_Dim_y(testCase)
	%Description:
	%	Verifies if Dim_y correctly identifies that a one-dimensional sytem's output is one-dimensional.

	%% Including Libraries
	include_relevant_libraries();

	%% Constants
	sys1 = get_1d_lcsas();

	%% Algorithm
	assert(	sys1.Dim_y() == 1 )

end

function test2_Dim_y(testCase)
	%Description:
	%	Verifies if Dim_y correctly identifies that the spacecraft system's output is two-dimensional.

	%% Including Libraries
	include_relevant_libraries();

	%% Constants
	sys2 = get_spacecraft_lcsas1();

	%% Algorithm
	assert(	sys2.Dim_y() == 2 )

end

function test3_Dim_y(testCase)
	%Description:
	%	Verifies if Dim_y correctly identifies that the 2 x 2 consensus system's output is 2*2*2-dimensional.

	%% Including Libraries
	include_relevant_libraries();

	%% Constants
	sys3 = get_consensus_dyn(2,2,0.1);

	%% Algorithm
	assert(	sys3.Dim_y() == 4*2 )

end

function test1_Dim_w(testCase)
	%Description:
	%	Verifies if Dim_w correctly identifies that a one-dimensional sytem's process noise is one-dimensional.

	%% Including Libraries
	include_relevant_libraries();

	%% Constants
	sys1 = get_1d_lcsas();

	%% Algorithm
	assert(	sys1.Dim_w() == 1 )

end

function test2_Dim_w(testCase)
	%Description:
	%	Verifies if Dim_w correctly identifies that the spacecraft system's process noise is two-dimensional.

	%% Including Libraries
	include_relevant_libraries();

	%% Constants
	sys2 = get_spacecraft_lcsas1();

	%% Algorithm
	assert(	sys2.Dim_w() == 2 )

end

function test3_Dim_w(testCase)
	%Description:
	%	Verifies if Dim_w correctly identifies that the 2 x 2 consensus system's process noise is two-dimensional.

	%% Including Libraries
	include_relevant_libraries();

	%% Constants
	sys3 = get_consensus_dyn(2,2,0.1);

	%% Algorithm
	assert(	sys3.Dim_w() == 2 )

end

function test1_Dim_v(testCase)
	%Description:
	%	Verifies if Dim_v correctly identifies that a one-dimensional sytem's output noise is one-dimensional.

	%% Including Libraries
	include_relevant_libraries();

	%% Constants
	sys1 = get_1d_lcsas();

	%% Algorithm
	assert(	sys1.Dim_v() == 1 )

end

function test2_Dim_v(testCase)
	%Description:
	%	Verifies if Dim_v correctly identifies that the spacecraft system's output noise is two-dimensional.

	%% Including Libraries
	include_relevant_libraries();

	%% Constants
	sys2 = get_spacecraft_lcsas1();

	%% Algorithm
	assert(	sys2.Dim_v() == sys2.Dim_y() )

end

function test3_Dim_v(testCase)
	%Description:
	%	Verifies if Dim_v correctly identifies that the 2 x 2 consensus system's output noise is 2*2*2-dimensional (same as the state, output).

	%% Including Libraries
	include_relevant_libraries();

	%% Constants
	sys3 = get_consensus_dyn(2,2,0.1);

	%% Algorithm
	assert(	sys3.Dim_v() == sys3.Dim_x() )

end