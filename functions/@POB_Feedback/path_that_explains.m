function path_index = path_that_explains( varargin )
	%path_that_explains
	%Description:
	%	This function finds which path from fb.LSequenceConsistencySets explains the observed external behavior.
	%
	%Usage:
	%	tf = fb.external_behavior_explained_by( external_behavior_in  )
	%	tf = fb.external_behavior_explained_by( external_behavior_in , 'verbosity' , 1 )

	%%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing  %%
	%%%%%%%%%%%%%%%%%%%%%%%

	pob_feedback = varargin{1};
	external_behavior_in = varargin{2};

	LSequenceConsistencySets = pob_feedback.LSequenceConsistencySets;

	argidx = 3;

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

	

	sequence_cset_contains_eb = [];
	for sequence_index = 1:length(LSequenceConsistencySets)
		sequence_cset_contains_eb = [ sequence_cset_contains_eb ; LSequenceConsistencySets(sequence_index).contains( external_behavior_in ) ];
	end

	%Return answer
	tf = all(  ); %If the problem is feasible, then the external behavior is explainable by the internal behavior set.

end