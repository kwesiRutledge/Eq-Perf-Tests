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
