%test_LCSAS.m
%Description:
%	Tests the function LCSAS.

function tests = test_LCSAS
	%disp(localfunctions)
	tests = functiontests(localfunctions);

function include_relevant_libraries()
	%Description:
	%	Attempts to add the relevant libraries/toolboxes to the path.

	%% Constants

	is_on_personal_mac = strcmp(getenv('USER'),'kwesirutledge');
	is_on_great_lakes = (strcmp(getenv('USER'),'krutledg') && isunix) ;

	addpath('../')

	%% Algorithm
	if is_on_personal_mac
		include_fcns2({'mosek','gurobi','tbxmanager'},'PathToDirectoryWithToolboxes','../../../')
	elseif is_on_great_lakes
		include_fcns2({'tbxmanager','YALMIP'},'PathToDirectoryWithToolboxes','../../../')
	else
		warning('Warning: There is no guaranteed performance for systems that are not on Kwesi''s computer.')
	end

