function [ constraint , new_variable ] = CreateNonemptyConstraint( ibs )
	%Description:
	%	Creates a YALMIP constraint that is satisfied when the InternalBehaviorSet ibs
	%	is not empty.

	%% Constants

	ibs_settings = ibs.ibs_settings;

	%% Algorithm

	switch ibs_settings.OpenLoopOrClosedLoop
	case 'Open'
		% Create Optimization Variable
		new_variable = sdpvar(ibs.Dim,1,'full');

		% Create Constraint
		constraint = [ ibs.A * new_variable <= ibs.b ] + [ ibs.Ae * new_variable == ibs.be  ];

	case 'Closed'
		% Create Closed Loop Matrices
		


	otherwise
		error(['Unexpected value of ibs_settings.OpenLoopOrClosedLoop: ' ibs_settings.OpenLoopOrClosedLoop ])
	end


end