%mpcMatrixTest.m

%% Constants

A = [1 1; 0 1];
B = [0;1];
C = [1,0];
f = [1;2];

try
	temp_sys1 = Aff_Dyn(A,B,f,C);
	temp_sys2 = Aff_Dyn(2,3,1,0.5);
catch
	error('You need to include the entire ''functions'' directory before running this test. \n Please update your path and try again.')
end

%% Test 1: Generating H from an Aff_Dyn Object

%Constants
A = [1 1; 0 1];
B = [0;1];
C = [1,0];
f = [1;2];

temp_sys1 = Aff_Dyn(A,B,f,C);

dim1.n_x = size(temp_sys1.A,1);
dim1.n_w = size(temp_sys1.B_w,2);

T = 3;
H = calc_w_effect_mat(temp_sys1,ones(T,1));
assert(isequal(size(H),[dim1.n_x*(T+1),dim1.n_x*T ]))

%% Test 2: Generating H from a struct

n = 3;
A = magic(n);

T = 3;
H = calc_w_effect_mat(A,T);
assert(isequal(size(H),[n*(T+1),n*T]));

%% Test 3: Correct Construction of H for Switched System

%Constants
A1 = [1 1; 0 1];
A2 = [1, 0; 1, 1];
A3 = zeros(2);
B = [0;1];
C = [1,0];
f = [1;2];

temp_sys_arr = [ Aff_Dyn(A1,B,f,C), Aff_Dyn(A2,B,f,C), Aff_Dyn(A3,B,f,C) ];
L = [1,2,3,1];
H = calc_w_effect_mat(temp_sys_arr,L);
assert(isequal(H(end-1:end,1:2),zeros(2)))
assert(isequal(H(5:6,1:2),A2))
assert(isequal(H(end-1:end,end-3:end-2),A1))

%% Test 4: Generating S from a struct

A = [1,1;0,1];
B = 0.5*ones(2,1);
T = 3;

S = calc_u_effect_mat(A,B,T);

%If constructed properly, the following submatrix should be within S.
assert(isequal(S(end-1:end,1),A^2*B)) 

%% Test 5a: Generating S from a Switched System (Part 1)

clear A1 A2 A3

A = [1 1; 0 1];
B1 = [0;1];
B2 = [1;0];
B3 = zeros(2,1);
C = [1,0];
f = [1;2];

temp_sys_arr = [ Aff_Dyn(A,B1,f,C), Aff_Dyn(A,B2,f,C), Aff_Dyn(A,B3,f,C) ];
L = [1,2,3,1];
S = calc_u_effect_mat(temp_sys_arr,L);

assert(isequal(S(end-1:end,end),B1))

%% Test 5b: Generating S from a Switched System (Part 2)
A = [1 1; 0 1];
B1 = [0;1];
B2 = [1;0];
B3 = zeros(2,1);
C = [1,0];
f = [1;2];

temp_sys_arr = [ Aff_Dyn(A,B1,f,C), Aff_Dyn(A,B2,f,C), Aff_Dyn(A,B3,f,C) ];
L = [1,2,3,1];
S = calc_u_effect_mat(temp_sys_arr,L);

%If constructed properly, the following submatrix should be within S.
assert(isequal(S(end-3:end-2,end-2),A*B2))

%% Test 6: Generating J from a Struct/Matrices

A = [1 1; 0 1];
x0 = [0;1];
T = 3;

J = calc_x0_mat(A,x0,T);
assert(isequal(J(end-1:end,1),(A^T)*x0))

%% Test 7a: Generating J from a Switched System (Part a)

%Constants
A1 = [1 1; 0 1];
A2 = [1, 0; 1, 1];
A3 = zeros(2);
L = [1,2,3];

B = [0;1];
C = [1,0];
f = [1;2];

temp_sys_arr = [ Aff_Dyn(A1,B,f,C), Aff_Dyn(A2,B,f,C), Aff_Dyn(A3,B,f,C) ];
L = [1,2,3];
J = calc_x0_mat(temp_sys_arr,L);

assert(isequal(J(end-1:end,1),zeros(2,1)))

%% Test 7b: Generating J from a Switched System (Part b)

%Constants
A1 = [1 1; 0 1];
A2 = [1, 0; 1, 1];
A3 = zeros(2);
L = [1,2,3];

B = [0;1];
C = [1,0];
f = [1;2];

temp_sys_arr = [ Aff_Dyn(A1,B,f,C), Aff_Dyn(A2,B,f,C), Aff_Dyn(A3,B,f,C) ];
L = [1,2,3];
temp_sys_arr(L(1)).x0 = [0;1];
J = calc_x0_mat(temp_sys_arr,L);

assert(isequal(J(end-3:end-2,1),A2*A1*[0;1]))

%% Test 8a: Full get_mpc_matrices test (Part a)

%Constants
A1 = [1 1; 0 1];
A2 = [1, 0; 1, 1];
A3 = zeros(2);
L = [1,2,3];

B = [0;1];
C = [1,0];
f = [1;2];

temp_sys_arr = [ Aff_Dyn(A1,B,f,C), Aff_Dyn(A2,B,f,C), Aff_Dyn(A3,B,f,C) ];
L = [1,2,3];
T = length(L);
[H,S,C_bar,J,f_bar] = get_mpc_matrices(temp_sys_arr,'time_horizon',T);

assert(isequal(size(f_bar),[2*T,1]))

%% Test 8b: Full get_mpc_matrices test (Part b)

%Constants
A1 = [1 1; 0 1];
A2 = [1, 0; 1, 1];
A3 = zeros(2);
L = [1,2,3];

B = [0;1];
C = [1,0];
f = [1;2];

temp_sys_arr = [ Aff_Dyn(A1,B,f,C), Aff_Dyn(A2,B,f,C), Aff_Dyn(A3,B,f,C) ];
L = [1,2,3];
T = length(L);
[H,S,C_bar,J,f_bar] = get_mpc_matrices(temp_sys_arr,'word',L);

assert(isequal(J(end-1:end,:),zeros(2)))

%% Test 8c: Full get_mpc_matrices test (Part c)

%Constants
A1 = [1 1; 0 1];
A2 = [1, 0; 1, 1];
A3 = zeros(2);
L = [1,2,3];

B = [0;1];
C = [1,0];
f = [1;2];

temp_sys_arr = [ Aff_Dyn(A1,B,f,C), Aff_Dyn(A2,B,f,C), Aff_Dyn(A3,B,f,C) ];
L = [1,2,3];
temp_sys_arr(L(1)).x0 = [0;1];
[H,S,C_bar,J,f_bar] = get_mpc_matrices(temp_sys_arr,'word',L);

assert(isequal(J(end-3:end-2,:),A2*A1));