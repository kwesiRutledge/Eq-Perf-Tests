function [results] = observer_tests( varargin )
	% 	observer_tests
	%		How to run the observer_comparisonX.m scriipts programmatically
	%	Usage:
	%		ot = observer_tests(80)
	%		ot = observer_tests(90,{'false,5'})
	%		ot = observer_Tests(90,{'false,5'},{'mosek','gurobi','tbxmanager'})

	%% Input Processing

	test_nums = varargin{1};

	if (nargin < 2)
		test_inputs = cell(length(varargin{1}),1);
	elseif (length(varargin{2})==0)
		test_inputs = cell(length(varargin{1}),1);
	else
		test_inputs = varargin{2};
	end

	% Determine which tools to include in path
	if nargin >= 3
		tools_to_include = varargin{3};
	else
		is_on_personal_mac = strcmp(getenv('USER'),'kwesirutledge');
		is_on_great_lakes = (strcmp(getenv('USER'),'krutledg') && isunix) ;

		%include_fcns('tbxmanager','YALMIP')
		if is_on_personal_mac
			tools_to_include = {'mosek','gurobi','tbxmanager'};
		elseif is_on_great_lakes
			tools_to_include = {'tbxmanager','YALMIP'};
		end
	end

	%% Include toolboxes.
	for tool_index = 1:length(tools_to_include)
		include_fcns(tools_to_include{tool_index});
	end

	%% Algorithm

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

	% Run tests provided in test_nums (can be array or integer)
    results = cell(length(test_nums),1);
	for k = 1 : length(test_nums)
		results{k} = eval([base_name num2str(test_nums(k)) '(' test_inputs{k} ')' ]);
	end

end