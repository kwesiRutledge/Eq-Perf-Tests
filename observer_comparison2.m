% observer_comparison2.m
%		This scriipt will compare the observers found using 3 different methods:
%			- Triangle Inequality
%			- Skaf and Boyd's Method
%			- Yong's Method
%		where the objective is to achieve Equalized Performance

clear all;
close all;
clc;

%%%%%%%%%%%%%%%
%% Constants %%
%%%%%%%%%%%%%%%

pl = 1;		%Performance Level
verbosity = 0;

%Import the Proper ACC Matrices

load('data/system_examples/acc_p.mat');

%Load appropriate functions
addpath('./functions/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First Comparison, Next Step Error %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% e0 		= sdpvar(size(acc.A,1),1,'full');
% delta  	= sdpvar(size(acc.E,2),1,'full');
% mu		= sdpvar(size(acc.C,1),1,'full');

% Triangle Inequality

L 		= sdpvar(size(acc.A,1),size(acc.C,1),'full');

alpha_0 = sdpvar(1,1,'full');

%Create optimization objective:

e1_t = norm(acc.A-L*acc.C,Inf) + norm(L,Inf) * (acc.m/pl) + norm(acc.E,Inf)*(acc.d/pl);

%Optimize

sol_ti = optimize([],e1_t,sdpsettings('verbose',verbosity));

if sol_ti.problem == 0
	disp('Optimization solved successfully.')
else
	error('Optimization error.')
end

% Skaf and Boyd Optimization
acc_e = acc;
acc_e.B = eye(size(acc.A,1));	%For error system, the input we are interested in
								% isn't modified by a B matrix
acc_e.x0 = Inf;		%The initial condition needs to be set due to my rigid,
					% old way of doing things. It will not be used when pl is given.

sol_skaf = generate_skaf_controller(acc_e,1,verbosity,'PL',pl);

%Create robust constraints on the Initial Condition, as well as the noise

% Yong Optimization

sol_yong = generate_yong_controller(acc_e,1,verbosity,'PL',pl);

save('results/oc_exp1.mat','sol_ti','sol_skaf','sol_yong');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Second Comparison, k^th Step Error %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear L alpha_0 e1_t

%Constants
t_lim = 10;

%Looping
for t =  1 : 10

	%Calculate Skaf and Yong Solutions at this time step
	sol_skaf = generate_skaf_controller(acc_e,t,verbosity,'PL',pl);
	sol_yong = generate_yong_controller(acc_e,t,verbosity,'PL',pl);

	%Save data to vector
	skaf_opt(t) = sol_skaf.opt_obj;
	yong_opt(t) = sol_yong.opt_obj;

end

figure;
hold on;
plot(skaf_opt)
plot(yong_opt)

legend('Skaf','Yong')