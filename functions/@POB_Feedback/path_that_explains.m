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

	% Output Measurements
	y_vec = reshape(y_mat,prod(size(y_mat)),1);
	u_vec = reshape(contr.u_hist,prod(size(contr.u_hist)),1);

	%Get the current time and find the nodes that correspond to this time.
	t = size(y_mat,2) - 1;

	sequence_cset_contains_eb = [];
	for sequence_index = 1:length(LSequenceConsistencySets)
		sequence_cset_contains_eb = [ sequence_cset_contains_eb ; LSequenceConsistencySets(sequence_index).contains( external_behavior_in ) ];
	end

	% Find the sequence that contained the external behavior and also had maximum cardinality.
	
	tf = all(  ); %If the problem is feasible, then the external behavior is explainable by the internal behavior set.

end