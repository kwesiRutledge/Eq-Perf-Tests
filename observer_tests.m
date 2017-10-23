function [results] = observer_tests( varargin )
	% 	observer_tests
	%		How to run the observer_comparisonX.m scriipts programmatically
	%

	if ~any(strcmp(path,'./functions'))
		addpath('./functions/')
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
		results(k) = eval([base_name num2str(test_nums(k)) '(' test_inputs{k} ')' ])
	end

end