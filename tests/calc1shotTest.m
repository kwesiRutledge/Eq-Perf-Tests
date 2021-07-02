% calc_1shot_isTests.m

%% Test 1

A_tilde = [1 1; 0 1];
B = [0.5;1];
K = [-0.0796, -0.4068];

A = A_tilde + B*K; % Cast as an autonomous linear system

F = [1 0; 0 1; -1 0; 0 -1];
g = 0.1*ones(4, 1);

P_w = Polyhedron('A',F,'b',g);

[~,info] = calc_1shot_is(A,P_w,20);

assert(info.exit_flag == 1)