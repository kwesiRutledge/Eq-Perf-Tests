% fixing_skaf.m

clear all;
close all;
clc;

%% Constants

x0 = randn(10,1);

%% Optimization Variables

Q = sdpvar(10,10,'full');
x = sdpvar(10,1,'full');

%% Constraints

%Lower Diagonal
l_diag_constr = [];
for row_num = 1 : size(Q,1)-1
	l_diag_constr = l_diag_constr + [ Q(row_num,[row_num+1:end]) == 0 ];
end

%Constraints on entry magnitude
entry_constrs = [];
for row_ind = 1:size(Q,1)
	for col_ind = 1 : size(Q,2)
		entry_constrs = entry_constrs + [ -1 <= Q( row_ind , col_ind ) <= 1 ];
	end
end

%Constraints on robust optimization variables.
robust_constrs = [];
robust_constrs = [ -3 <= x <= 3, uncertain(x) ];



%% Optimize


sol = optimize(l_diag_constr + entry_constrs , norm(Q*x + (eye(size(Q)) + Q)*x0,Inf) , sdpsettings('verbose',1) );

if sol.problem == 0
	disp(['YALMIP Optimization Solved'])
else
	%error(['YALMIP Robust Optimization #' num2str(T) ' NOT Solved.'])
	disp(['YALMIP Optimization NOT Solved. (' yalmiperror(sol.problem) ').'])
end

disp([ 'objective value: ' num2str(value(norm(Q,Inf)) ) ])