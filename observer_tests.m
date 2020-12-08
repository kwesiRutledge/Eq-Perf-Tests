function [results] = observer_tests( varargin )
	% 	observer_tests
	%		How to run the observer_comparisonX.m scriipts programmatically
	%

	% Add experiments and functions to the path
	is_on_personal_mac = strcmp(getenv('USER'),'kwesirutledge');
	is_on_great_lakes = (strcmp(getenv('USER'),'krutledg') && isunix) ;

	%include_fcns('tbxmanager','YALMIP')
	if is_on_personal_mac
		include_fcns('mosek','gurobi','tbxmanager')
	elseif is_on_great_lakes
		include_fcns('tbxmanager','YALMIP')
	end

	if isempty(strfind(path,'./functions/'))
		addpath('./functions/')
		addpath('./functions/systems/')
	end

	try
		empty_function_experiments()
	catch
		addpath('./experiments/')
	end

	%% Constants
	base_name = 'observer_comparison';

	test_nums = varargin{1};
	if nargin < 2
		test_inputs = cell(length(varargin{1}),1);
	else
		test_inputs = varargin{2};
	end

	%%Run tests from 
    results = cell(length(test_nums),1);
	for k = 1 : length(test_nums)
		results{k} = eval([base_name num2str(test_nums(k)) '(' test_inputs{k} ')' ]);
	end

end