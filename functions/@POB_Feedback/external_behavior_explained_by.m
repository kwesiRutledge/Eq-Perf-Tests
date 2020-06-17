function tf = external_behavior_explained_by( varargin )
	%external_behavior_explained_by
	%Description:
	%	This function evaluates whether or not an external behavior (external_behavior_in) can be generated by
	%	an internal behavior from internal_behavior_set_in.
	%	POB_Feedback.external_behavior_explained_by()
	%
	%Usage:
	%	tf = fb.external_behavior_explained_by( external_behavior_in , internal_behavior_set_in )
	%	tf = fb.external_behavior_explained_by( external_behavior_in , internal_behavior_set_in , 'verbosity' , 1 )

	%%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing  %%
	%%%%%%%%%%%%%%%%%%%%%%%

	pob_feedback = varargin{1};
	external_behavior_in = varargin{2};
	internal_behavior_set_in = varargin{3};

	argidx = 4;

	while argidx <= nargin
		switch varargin{argidx}
			case 'verbosity'
				verbosity = varargin{argidx+1};
				argidx = argidx + 2;
			otherwise
				error('The input to external_behavior_explained_by() is not available.')
		end
	end

	if ~exist('verbosity')
		verbosity = 0;
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	ib_dim = internal_behavior_set_in.Dim;
	eb_dim = length(external_behavior_in);

	options = sdpsettings('verbose',verbosity);

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	% Construct Optimization Variables
	internal_behavior = sdpvar(ib_dim,1,'full');


	% Construct Constraints
	ib_match_constraint = [ [eye(eb_dim) , zeros(eb_dim,ib_dim-eb_dim)]*internal_behavior == external_behavior_in ];

	valid_ib_constraint = 	[ internal_behavior_set_in.A * internal_behavior <= internal_behavior_set_in.b ] + ...
							[ internal_behavior_set_in.Ae * internal_behavior == internal_behavior_set_in.be ];

	%Run optimization
	diagnostics = optimize(	ib_match_constraint+valid_ib_constraint, ...
							[], options);

	%Return answer
	tf = (diagnostics.problem == 0); %If the problem is feasible, then the external behavior is explainable by the internal behavior set.

end