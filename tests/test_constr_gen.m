%test_constr_gen.m
%Description:
%	This script is meant to test various functions for the constr_gen.m class in this repository.

function tests = test_get_1d_lcsas
	%disp(localfunctions)
	tests = functiontests(localfunctions);
end

function include_relevant_libraries()
	%Description:
	%	Attempts to add the relevant libraries/toolboxes to the path.

	%% Constants

	is_on_personal_mac = strcmp(getenv('USER'),'kwesirutledge');
	is_on_great_lakes = (strcmp(getenv('USER'),'krutledg') && isunix) ;

	%% Algorithm

	%Include Yalmip
	%include_fcns2({'mosek','gurobi','MPT3'})

	%Add Local Functions to Path
	addpath(genpath('../functions'));
end

function test_binary_inclusion_constraint1(testCase)
	%Description:
	%	This function is meant to create an constraint and the accompanying optimization variables such that
	%		"there exists an x 		such that constraint  Ax <= b        		hold, if and only if the binary variable binvar = 1" and
	%		"there exists a y >= 0 	such that constraints b^Ty < 0 and A^Ty = 0 hold, if and only if the binary variable binvar = 0"

	%% Include Relevant Libraries
	% include_relevant_libraries();

	%% Constants

	n = 2;
	P1 = Polyhedron('lb',-ones(1,n),'ub',ones(1,n));
	P2 = Polyhedron('A',[-1,0;0,-1;-1,1],'b',[0;0;1]);

	eps0 = 10^(-3);
	ops0 = sdpsettings('verbose',0);

	%% Algorithm 

	A = [P1.A;P2.A];
	b = [P1.b;P2.b];

	mA = size(A,2);
	nA = size(A,1);

	iv1 = binvar(1,1);
	x = sdpvar(mA,1,'full');
	y = sdpvar(nA,1,'full');

	pos_constr = [ y >= 0 ];

	constr1 = [ iv1 .* A * x + (1-iv1).*(b'*y + eps0) <= iv1 .* b + (1 - iv1) .* 0 ];
	constr2 = [  (1-iv1)*A'*y == 0 ];

	% Solve the optimization problem

	opt_out = optimize( pos_constr + constr1 + constr2 , [] , ops0 );

	disp('Results:')
	disp(['- iv1 = ' num2str(value(iv1)) ])
	disp(['- x'' = ' num2str(value(x')) ])
	disp(['- y'' = ' num2str(value(y')) ])

	assert( (opt_out.problem == 0) && (value(iv1) == 1) )

end

function test_binary_inclusion_constraint2(testCase)
	%Description:
	%	This function is meant to create an constraint and the accompanying optimization variables such that
	%		"there exists an x 		such that constraint  Ax <= b        		hold, if and only if the binary variable binvar = 1" and
	%		"there exists a y >= 0 	such that constraints b^Ty < 0 and A^Ty = 0 hold, if and only if the binary variable binvar = 0"
	%	In this example, there should not exist an x (i.e. b = 0).

	%% Include Relevant Libraries
	% include_relevant_libraries();

	%% Constants

	n = 2;
	P1 = Polyhedron('lb',-ones(1,n),'ub',ones(1,n));
	P2 = Polyhedron('lb',[1.5,1.5],'ub',[2,2]);

	eps0 = 10^(-3);
	ops0 = sdpsettings('verbose',0);

	%% Algorithm 

	A = [P1.A;P2.A];
	b = [P1.b;P2.b];

	mA = size(A,2);
	nA = size(A,1);

	iv1 = binvar(1,1);
	x = sdpvar(mA,1,'full');
	y = sdpvar(nA,1,'full');

	pos_constr = [ y >= 0 ];

	constr1 = [ iv1 .* A * x + (1-iv1).*(b'*y + eps0) <= iv1 .* b + (1 - iv1) .* 0 ];
	constr2 = [  (1-iv1)*A'*y == 0 ];

	% Solve the optimization problem

	opt_out = optimize( pos_constr + constr1 + constr2 , [] , ops0 );

	disp('Results:')
	disp(['- iv1 = ' num2str(value(iv1)) ])
	disp(['- x'' = ' num2str(value(x')) ])
	disp(['- y'' = ' num2str(value(y')) ])

	assert( (opt_out.problem == 0) && (value(iv1) == 0) )

end

function test_binary_inclusion_constraint3(testCase)
	%Description:
	%	This function is meant to create an constraint and the accompanying optimization variables such that
	%		"there exists an x 		such that constraint  Ax <= b        		hold, if and only if the binary variable binvar = 1" and
	%		"there exists a y >= 0 	such that constraints b^Ty < 0 and A^Ty = 0 hold, if and only if the binary variable binvar = 0"
	%
	%

	%% Include Relevant Libraries
	include_relevant_libraries();

	%% Constants

	n = 2;
	P1 = Polyhedron('lb',-ones(1,n),'ub',ones(1,n));
	P2 = Polyhedron('A',[-1,0;0,-1;-1,1],'b',[0;0;1]);

	ops0 = sdpsettings('verbose',0);

	cg = constr_gen();

	%% Algorithm 

	A = [P1.A;P2.A];
	b = [P1.b;P2.b];

	[ constrs , iv1 , x , y ] = cg.binary_inclusion_constraint(A,b);

	% Solve the optimization problem

	opt_out = optimize( constrs , [] , ops0 );

	disp('Results:')
	disp(['- iv1 = ' num2str(value(iv1)) ])
	disp(['- x'' = ' num2str(value(x')) ])
	disp(['- y'' = ' num2str(value(y')) ])

	assert( (opt_out.problem == 0) && (value(iv1) == 1) )

end