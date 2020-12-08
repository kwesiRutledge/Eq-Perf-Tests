function [constraints,binary_var,x,y] = binary_inclusion_constraint( obj , A , b )
	%Description:
	%	This creates constraints, a binary variable and continuous variables such that use Farkas Lemma such that 

	%% Input Processing %%

	%% Constants %%

	mA = size(A,2);
	nA = size(A,1);

	eps0 = 10^(-3);

	%% Create Optimization Variables %%

	iv1 = binvar(1,1);
	x = sdpvar(mA,1,'full');
	y = sdpvar(nA,1,'full');

	%% Create Constraints %%

	pos_constr = [ y >= 0 ];

	constr1 = [ iv1 .* A * x + (1-iv1).*(b'*y + eps0) <= iv1 .* b + (1 - iv1) .* 0 ];
	constr2 = [  (1-iv1)*A'*y == 0 ];

	%% Create Outputs %%

	binary_var = iv1;
	constraints = pos_constr + constr1 + constr2;


end