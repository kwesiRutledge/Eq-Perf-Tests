function [ o1, constr ] = create_zonotope_inclusion_constraint( )

	%Create sdpvar
	o1 = sdpvar(1,1,'full');

	constr = [o1 <= 2];

end