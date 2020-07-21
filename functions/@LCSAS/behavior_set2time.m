function [t_out] = behavior_set2time( varargin )
	%behavior_set2time.m
	%Description:
	%
	%
	%Usage:
	%	t_out = lcsas0.behavior_set2time( set_in , 'ConsistencySet' )
	%	t_out = lcsas0.behavior_set2time( set_in , 'InternalBehaviorSet_1word' )
	%	t_out = lcsas0.behavior_set2time( set_in , 'ConsistencySet' , 'ConsistencySetVersion' , 1 )
	%	t_out = lcsas0.behavior_set2time( set_in , 'ConsistencySet' , 'FeedbackMethod' , 'state' )
	%

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	lcsas_in = varargin{1};
	set_in = varargin{2};
	set_type = varargin{3};

	arg_idx = 4;
	while arg_idx <= nargin
		switch varargin{arg_idx}
			case 'ConsistencySetVersion'
				ConsistencySetVersion = varargin{arg_idx+1};
				arg_idx = arg_idx + 2;
			case 'FeedbackMethod'
				FeedbackMethod = varargin{arg_idx+1};
				arg_idx = arg_idx + 2;
			otherwise
				error('Unrecognized input to this function.')
		end
	end

	if ~exist('ConsistencySetVersion')
		ConsistencySetVersion = 1;
	end

	if ~exist('FeedbackMethod')
		FeedbackMethod = 'output';
	end

	if ~any(strcmp( FeedbackMethod , {'state','output'} ) )
		error(['Unrecognized choice for FeedbackMethod: ' FeedbackMethod] )
	end

	set_type_strings = {'ConsistencySet','InternalBehaviorSet_1word'};
	if ~any(strcmp( set_type , set_type_strings ))
		error(['Unrecognized choice for set_type: ' set_type ])
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	n_x = size(lcsas_in.Dyn(1).A,1);
	n_u = size(lcsas_in.Dyn(1).B,2);
	n_w = size(lcsas_in.Dyn(1).B_w,2);
	n_y = size(lcsas_in.Dyn(1).C,1);
	n_v = size(lcsas_in.Dyn(1).C_v,2);

	dim = set_in.Dim;

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	ConsistencySetVersion;
	
	switch set_type
		case 'ConsistencySet'
			if strcmp(FeedbackMethod,'output')
				t_out = (dim-n_y)/(n_y+n_u);
			elseif strcmp(FeedbackMethod,'state')
				t_out = (dim - n_x)/(n_x+n_u);
				warning('Untested use of behavior_set2time!')
			end
		case 'InternalBehaviorSet_1word'
			%This must be handled depending on the version of the consistency set
			switch ConsistencySetVersion
				case 1
					if strcmp(FeedbackMethod,'output')
						t_out = (dim-n_y-n_v-2*n_x)/(n_y+n_u+n_w+n_v+n_x);
					elseif strcmp(FeedbackMethod,'state')
						t_out = (dim - 2*n_x)/(n_x+n_u+n_w);
					end
				case 2
					if strcmp(FeedbackMethod,'output')
						t_out = (dim - n_y - n_x)/(n_y+n_u+n_w);
					elseif strcmp(FeedbackMethod,'state')
						error('This part is not implemented!')
					end
				otherwise
					error(['The ConsistencySetVersion value ' num2str(ConsistencySetVersion) ' is not valid.' ])
			end
		otherwise
			error(['Unrecognized set_type (' set_type ').'])
	end
		


	if (round(t_out) ~= t_out)
		error(['Your choice of flags may be incorrect. t is not an integer (' num2str(t_out) ').'])
	end

end