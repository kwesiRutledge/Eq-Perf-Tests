function [results] = observer_comparison34(varargin)
%observer_comparison34.m
%   This is an outgrowth of the results of experiment 33, which hinted that
%   when the system is "observable" and we have a long enough time horizon 
%   T, there exists a certain pair of parameters m1-bar and m3-bar such
%   that the problem (m1-bar,m2,m3-bar) is feasible AND for all m1 > m1-bar
%   the problem (m1,m2,m3-bar) is feasible.

%% Constants

L = [ 1, 1, 0 ];

eta_w = 0.1;
eta_v = 0.2;

M1 = 1;

expm_under_test = [2];

%% Creating a few interesting systems to test

A = [ 1 0 0;
      1 1 0;
      -1 5 -2];

[v1,d1] = eig(A)
  
B = [ 1; 1; 0 ];
C = [0 1 0]; % v1(:,[1 3])';

K_sub = place(A([1 2],[1 2]),B([1 2],:),[-0.5 0.5]);

ad1 = Aff_Dyn(A,B,zeros(size(A,1),1),C,...
                eta_w,eta_v);

%% Synthesizing Best Guarantees for the problem with Detectable System as M1 varies
if any(expm_under_test == 1)
    M3_vals = [];
    for M1_val = [0:0.1:2]
        [opt_out,estima] = ad1.free_rec_design_pb('Min_M3',M1_val,100,L,'verbosity',0);

        M3_vals = [M3_vals opt_out.M3];

        fprintf('.')
    end

    exp1.M1 = [0:0.1:2];
    exp1.M3 = M3_vals;

    figure;
    plot(exp1.M1([1:length(M3_vals)]),exp1.M3)
    xlabel('M1 Assumption Value')
    ylabel('Minimum M3 Value')
end

%% Synthesizing Optimal M3 Guarantee for Problem 2

%Experiment 2 Parameters
%   - Detectible System (A,C)
%   - Input Matrix B is identity

if any(expm_under_test == 2)

    ad2 = ad1;
    ad2.B = eye(3);

    M3_vals = [];
    for M1_val = [0:0.1:2]
        [opt_out,estima] = ad2.free_rec_design_pb('Min_M3',M1_val,100,L,'verbosity',0);

        M3_vals = [M3_vals opt_out.M3];

        fprintf('.')
    end
    
    exp1.M1 = [0:0.1:2];
    exp1.M3 = M3_vals;

    figure;
    plot(exp1.M1([1:length(M3_vals)]),exp1.M3)
    xlabel('M1 Assumption Value')
    ylabel('Minimum M3 Value')
    title('Experiment 2: Identity B Matrix')
    
end

%% Assigning Results

results.exp1 = exp1;

end