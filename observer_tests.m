function [results] = observer_tests( varargin )
	% 	observer_tests
	%		How to run the observer_comparisonX.m scriipts programmatically
	%

	% Add experiments and functions to the path
	%include_fcns('tbxmanager','YALMIP')
	include_fcns('mosek','gurobi','MPT3')

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
		test_inputs = cell(size(varargin{1}));
	else
		test_inputs = varargin{2};
	end

	%%Run tests from 
	for k = 1 : length(test_nums)
		results{k} = eval([base_name num2str(test_nums(k)) '(' test_inputs{k} ')' ]);
	end

end