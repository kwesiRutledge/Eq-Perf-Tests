function [ constraint , new_variable ] = CreateNonemptyConstraint( ibs )
	%Description:
	%	Creates a YALMIP constraint that is satisfied when the InternalBehaviorSet ibs
	%	is not empty.

	%% Constants

	ibs_settings = ibs.ibs_settings;

	%% Algorithm

	% Create Optimization Variable
	new_variable = sdpvar(ibs.Dim,1,'full');

	% Create Constraint
	constraint = [ ibs.A * new_variable <= ibs.b ] + [ ibs.Ae * new_variable == ibs.be  ];


end