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

function test_create_mccormick_envelope1(testCase)
	%test_create_mccormick_envelope1
	%Description:
	%	Tests the error catching capabilities of the create_mccormick_envelope function.
	%	This test should recognize that the input of two matrix sdpvar objects is invalid/not supported.

	%% Constants 

	cg = constr_gen(0);

	%% Algorithm 
	x = sdpvar(2,2,'full');
	y = sdpvar(2,2,'symmetric');

	try
		cg.create_mccormick_envelope(x,y,[])
	catch e
		% disp(e.message)
		assert(strcmp(e.message,'It must be the case that x and y are of compatible dimensions. Make sure that x,y are both scalars, both vectors, or a matrix and a vector.'))
	end

end

function test_create_mccormick_envelope2(testCase)
	%test_create_mccormick_envelope2
	%Description:
	%	Tests the error catching capabilities of the create_mccormick_envelope function.
	%	This test should recognize that the length of two input vectors should be the same.

	%% Constants 

	cg = constr_gen(0);

	%% Algorithm 
	x = sdpvar(2,1,'full');
	y = sdpvar(3,1,'full');

	tempGridDef = struct('x_lb',0,'x_ub',1,'y_lb',0,'y_ub',1,'NumberOfRegionsInX',2,'NumberOfRegionsInY',3);

	try
		cg.create_mccormick_envelope(x,y,tempGridDef)
	catch e
		% disp(e.message)
		assert(strcmp(e.message,'Expected the lengths of x and y to be the same. Instead, length(x) = 2 and length(y) = 3.'))
	end

end

function test_create_mccormick_envelope3(testCase)
	%test_create_mccormick_envelope2
	%Description:
	%	Tests the error catching capabilities of the create_mccormick_envelope function.
	%	This test should recognize that the input struct does not have field NumberOfRegionsInX.

	%% Constants 

	cg = constr_gen(0);

	%% Algorithm 
	x = sdpvar(2,1,'full');
	y = sdpvar(3,1,'full');

	tempGridDef = struct('x_lb',0,'x_ub',1,'y_lb',0,'y_ub',1,'NumberOfRegionsInY',3);

	try
		cg.create_mccormick_envelope(x,y,tempGridDef)
	catch e
		%disp(e.message)
		assert(strcmp(e.message,[ 'The field name NumberOfRegionsInX is not part of the input grid_definition.' ]))
	end

end

function test_create_mccormick_envelope4(testCase)
	%Description:
	%	Tests to see if the expected answer of a known optimization problem matches the output
	%	of create_mccormick_envelope().

	%% Constants

	cg = constr_gen(0);
	eps0 = 10^(-1);

	num_ticks_in_K = 20;
	num_ticks_in_y = 20;

	%% Algorithm

	% Create Variables

	K = sdpvar(1,1,'full');
	y = sdpvar(1,1,'full');

	% z = sdpvar(1,1,'full');
	% mccormick_binvars = binvar(num_mccormick_envelopes,1,'full');

	% Create Constraints

	y_lb = 0;
	y_ub = 40;
	K_lb = eps0;
	K_ub = 1;

	range_constrs = [ K_lb <= K <= K_ub ] + [ y_lb <= y <= y_ub ];

	tempGridDef = struct('x_lb',K_lb,'x_ub',K_ub,'y_lb',y_lb,'y_ub',y_ub,'NumberOfRegionsInX',num_ticks_in_K,'NumberOfRegionsInY',num_ticks_in_y);

	[ z , mccormick_constraints , mccormick_binvars ] = cg.create_mccormick_envelope( K , y , tempGridDef , 'eta_z_bounds' , 10^5  );
	

	% ideal_constraint = [ K*y <= 1 ]

	bounding_constraint = [ z <= 1 ];

	%%%%%%%%%%%%%%%%
	%% Optimize ? %%
	%%%%%%%%%%%%%%%%

	optimization_constraints = range_constrs + mccormick_constraints + bounding_constraint;

	ops = sdpsettings('verbose',0,'debug',1);
	ops = sdpsettings(ops, 'solver','gurobi');
	optim0 = optimize(optimization_constraints,[-(y)],ops);

	assert( (value(y) > 10-eps0) && (value(y) < 10+eps0))

end

function test_create_mccormick_envelope5(testCase)
	%test_create_mccormick_envelope5
	%Description:
	%	Tests to see if the expected answer of a known optimization problem matches the output
	%	of create_mccormick_envelope().

	%% Constants

	cg = constr_gen(0);
	eps0 = 10^(-1);

	num_ticks_in_K = [10,10];
	num_ticks_in_w = [10,10];

	Pw1 = Polyhedron('lb',[0,2],'ub',[1,3]);
	Pw1.outerApprox;

	%% Algorithm

	% Create Variables

	K = sdpvar(1,2,'full');
	w = sdpvar(2,1,'full');

	% z = sdpvar(1,1,'full');
	% mccormick_binvars = binvar(num_mccormick_envelopes,1,'full');

	% Create Constraints

	w_lb = []; w_ub=[];
	for dim_idx = 1:Pw1.Dim 
		w_lb(dim_idx) = Pw1.Internal.lb(dim_idx);
		w_ub(dim_idx) = Pw1.Internal.ub(dim_idx);
	end

	K_lb = [ -10 ; -10 ];
	K_ub = [ 10 ; 10 ];

	range_constrs = [ K_lb <= K' <= K_ub ] + [ Pw1.A*w <= Pw1.b ];

	tempGridDef = struct('x_lb',K_lb,'x_ub',K_ub,'y_lb',w_lb,'y_ub',w_ub,'NumberOfRegionsInX',num_ticks_in_K,'NumberOfRegionsInY',num_ticks_in_w);

	[ z , mccormick_constraints , mccormick_binvars ] = cg.create_mccormick_envelope( K , w , tempGridDef , 'eta_z_bounds' , 10^5  );
	

	% ideal_constraint = [ K*y <= 1 ]

	bounding_constraint = [ z <= 1 ];

	%%%%%%%%%%%%%%%%
	%% Optimize ? %%
	%%%%%%%%%%%%%%%%

	optimization_constraints = range_constrs + mccormick_constraints + bounding_constraint;

	ops = sdpsettings('verbose',0,'debug',1);
	ops = sdpsettings(ops, 'solver','gurobi');
	optim0 = optimize(optimization_constraints,[-(w(2))],ops);

	assert( (value(w(2)) > 3-eps0) && (value(w(2)) < 3+eps0))

end

function test_get_H_polyt_inclusion_constr1(testCase)
	%Description:
	%	Verify if the symbolic version of the function works.

	H_x = [1,2;3,4];
	h_x = [5;6];
	H_y = [7,8;9,10];
	h_y = [11;12];

	cg = constr_gen(0);

	[ Lambda , eq_constr , ineq_constr ] = cg.get_H_polyt_inclusion_constr(H_x,h_x,H_y,h_y,'ReturnValueType','Symbolic')

	% Check All Equality Constraint Terms

	eq_row_index = 1;
	eq_col_index = 1;
	[c,t] = coeffs(eq_constr(eq_row_index,eq_col_index),Lambda);

	assert(all(c==[1,3,-7]) )

	eq_row_index = 2;
	eq_col_index = 1;
	[c,t] = coeffs(eq_constr(eq_row_index,eq_col_index),Lambda);

	assert(all(c==[1,3,-9]))

	eq_row_index = 1;
	eq_col_index = 2;
	[c,t] = coeffs(eq_constr(eq_row_index,eq_col_index),Lambda);

	assert(all(c==[2,4,-8]) )

	eq_row_index = 2;
	eq_col_index = 2;
	[c,t] = coeffs(eq_constr(eq_row_index,eq_col_index),Lambda);

	assert(all(c==[2,4,-10]))

	% Check all inequality constraints

	[c,t] = coeffs(ineq_constr(1,1),Lambda);

	assert(all(c==[5,6,-11]) )

	[c,t] = coeffs(ineq_constr(2,1),Lambda);

	assert(all(c==[5,6,-12]))


end