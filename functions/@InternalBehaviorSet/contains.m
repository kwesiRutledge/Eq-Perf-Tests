function tf = contains(ibs,internal_behavior)
	%Description:
	%	Determines whether internal_behavior is contained in (true) or not contained in (false) ibs.

	%% Constants %%

	ib_dim = ibs.Dim;
    eps0 = 10^(-5);

	%% Input Processing %%

	if ~isvector(internal_behavior)
		error(['internal_behavior input to contains() is not a vector!'])
	end

	if length(internal_behavior) ~= ib_dim
		error(['InternalBehaviorSet is of dimension ' num2str(ib_dim) ', but the input internal_behavior has dimension ' num2str(length(internal_behavior)) '.' ])
	end

	%% Algorithm %%

	tf = all( ibs.A * internal_behavior <= ibs.b ) ...
        && all( ibs.Ae * internal_behavior - ibs.be < eps0 ) && all( -(ibs.Ae * internal_behavior - ibs.be) < eps0 ) ; %Equality is approximately true

end