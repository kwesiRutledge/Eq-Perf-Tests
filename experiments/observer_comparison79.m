function [results] = observer_comparison79( varargin )
	%observer_comparison79.m
	%Description:
	%	Plotting data after running observer_comparison77.m for multiple conditions.
	%

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	% if nargin >= 1
	% 	experimental_post_flag = varargin{1};
	% end

	% if nargin >= 2
	% 	system_name = varargin{2};
	% end

	% if nargin >= 3
	% 	verbosity = varargin{3};
	% end

	% if nargin >= 3
	% 	c_sq.dim_x = varargin{2};
	% 	c_sq.dim_y = varargin{3};
	% end

	% if nargin >= 4
	% 	verbosity = varargin{4};
	% end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	twod_sys_res = [ mean([24,32.663,31.725]), mean([53,53.854,56.887]) , mean([57,26.622,22.968]) , NaN ];
	twoByOne_sys_res = [ 106.017, 833 , NaN , mean([958,1033])];
	compass_sys_res = [ 376 , 76.692 , NaN , 91.911 ];

	results.twod_sys_res = twod_sys_res;
	results.twoByOne_sys_res = twoByOne_sys_res;
	results.compass_sys_res = compass_sys_res;

	%% Plotting %%

	X = categorical({'2D Simple','2x1 Consensus','Compass Walker'});
	X = reordercats(X,{'2D Simple','2x1 Consensus','Compass Walker'});

	figure;
	bar( X , [ twod_sys_res([1,2])' , twoByOne_sys_res([1,2])' , compass_sys_res([1,2])' ])
	ylabel('Time (s)')

	twod_sys_res_mod1 = [ mean([24,32.663,31.725]), mean([53,53.854,56.887]) , mean([26.622,22.968]) , NaN ];
	results.twod_sys_res_modulo_outlier = twod_sys_res_mod1;

end