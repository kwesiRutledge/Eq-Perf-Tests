function [Inv,info] = calc_1shot_is(A,P_w,r)
%calc_1shot_is Approximates the invariant set using the 1 shot method
%developed in Paul Trodden
%   Detailed Description:
%       Uses Paul Trodden's 2016 work "A One-Step Approach to Computing a
%       Polytopic Robust Positively Invariant Set", we create a function
%       which should calculate the robust invariant set for discrete time
%       systems of the form:
%           x^+ = A x + w
%       where w is assumed to belong to a polyhedral set .
%
%   Inputs:
%       A   -   The state update matrix.
%       P_w -   The polyhedral set that the disturbance w belongs to.
%       r   -   The number of hyperplanes considered in Trodden's method.
%
%   Outputs:
%       Inv -   The invariant set for the system described by the inputs
%       and the above equation.


    %% Input Processing
    if size(A,1) ~= size(A,2)
        error('''A'' matrix is not square.')
    end
    
    %% Constants
    p = size(P_w.A,1);
    n = size(P_w.A,2);
    
    %% Algorithm
    % Predefine P as an r-sided regular polyhedron
    P = zeros(r, 2);
    for i = 1:r
        P(i,:) = [sin(2*pi*(i - 1)/r) cos(2*pi*(i - 1)/r)];
    end
    %P_poly = Polyhedron('A', P, 'b', ones(r, 1));
    %plot(P_poly)

    % Solve the linear program
    f = [ones(1, 2*r) zeros(1, 2*r*n)]; % Sum of c and d

    % Inequality matrix
    X = zeros(2*r + r^2 + r*p, 2*r + 2*n*r);

    % c_i - P_i*A*zeta_i <= 0
    X(1:r, 1:r) = eye(r);
    for i = 1:r
        X(i, ((i - 1)*n + 2*r + 1):(i*n + 2*r)) = -P(i, :)*A;
    end

    % -c - d + P*zeta_i <= 0
    for i = 1:r
        X(((i - 1)*r+(r+1)):((i*r) + r), 1:r) = -1*eye(r);
        X(((i - 1)*r+(r+1)):((i*r) + r), (r + 1):(2*r)) = -1*eye(r);
        X(((i - 1)*r+(r+1)):((i*r) + r), ((i - 1)*n + 2*r + 1):(i*n + 2*r)) = P;
    end

    % d_i - P_i*omega_i <= 0
    X((r^2 + r + 1):(r^2 + 2*r), (r+1):(2*r)) = eye(r);
    for i = 1:r
        X((r^2 + r + i), (2*r + n*r + (i - 1)*n + 1):(i*n + n*r + 2*r)) = -P(i, :);
    end

    % F*omega_i <= g
    for i = 1:r
        X((2*r + r^2 + (i - 1)*p + 1):(2*r + r^2 + i*p), (2*r + n*r + (i - 1)*n + 1): ...
            (2*r + n*r + i*n)) = P_w.A;
    end

    x = zeros(1, 2*r + r^2);
    for i = 1:r
        x = [x P_w.b'];
    end

    [sol,fval,exit_flag] = linprog(-f,X,x);
    q = sol(1:r) + sol((r+1):2*r);
    Inv = Polyhedron('A', P, 'b', q);
    
    info.fval = fval;
    info.exit_flag = exit_flag;

end

