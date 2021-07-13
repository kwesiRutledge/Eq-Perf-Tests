function tf = containsExternalBehavior( ibs , external_behavior_in )
	%Description:
	%
	%Usage:
	%	tf = ebs.contains( external_behavior )

	%% Constants

	ib_dim = ibs.Dim;

	%% Algorithm
	ebs = ibs.ToExternalBehaviorSet();

	tf = ebs.contains(external_behavior_in);

end