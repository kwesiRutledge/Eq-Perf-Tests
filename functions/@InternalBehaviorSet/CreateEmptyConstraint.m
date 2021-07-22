function [ constraint , new_variable ] = CreateEmptyConstraint(ibs)
	%Description:
	%	Creates a YALMIP constraint that is satisfied when the InternalBehaviorSet ibs
	%	is Empty.

	%% Constants	
	A = ibs.A;
	Ae = ibs.Ae;

	b = ibs.b;
	be = ibs.be;

	eps0 = 10^(-4);

	%% Create Constraints in YALMIP using Farkas Lemma
	H = [ A ; Ae ; -Ae ];
	h = [ b ; be ; -be ];

	n_H = size(H,1);
	new_variable = sdpvar(n_H,1,'full');

	constraint = [ H'*new_variable == 0] + [ h'*new_variable <= -eps0 ] + [new_variable >= 0];

end