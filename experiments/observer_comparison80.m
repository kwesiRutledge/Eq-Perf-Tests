function [results] = observer_comparison80( varargin )
	%observer_comparison80.m
	%Description:
	%	Playing around with the compass walker dynamics.
	%

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	test_name = 'observer_comparison80';

	disp(['Beginning ' test_name '.' ])
	disp(' ')

	disp('Defining Constants.')

	g = 10;

	if ~exist('a')
		a = 0.5; 	%meters
		b = a;

		l = a+b;
		m = 5; 		%kg
		m_H = 2*m; 	%kg
	end

	% theta_ns = x(1);
	% theta_s = x(2);
	% dot_theta_ns = x(3);
	% dot_theta_s = x(4);

	%Swing Phase Constants
	M = @(theta_ns,theta_s) [ m*(b^2) , -m*l*b*cos( theta_s - theta_ns ) ;
				-m*l*b*cos(theta_s - theta_ns) , (m_H+m)*l^2+m*a^2 ];
	
	N = @(theta_ns,theta_s,dot_theta_ns,dot_theta_s) ...
			[ 	0 , m*l*b*dot_theta_s*sin(theta_s - theta_ns);
				-m*l*b*dot_theta_ns*sin(theta_s - theta_ns) , 0 ];

	G = @(theta_ns,theta_s) [ m*b*sin(theta_ns) ; -(m_H*l+m*a+m*l)*sin(theta_s) ] * g;

	S = [	1,-1,0;
			0,1,1];

	%Transition Equation Constants
	if ~exist('alpha0')
		alpha0 = deg2rad(30);
	end

	P_alpha = [ (m_H*l^2 + 2*m*l^2)*cos(2*alpha0) - m*a*b - 2*m*b*l*cos(2*alpha0) , -m*a*b ;
				-m*a*b , 0 ];

	Q_alpha = [ m*b^2 - m*b*l*cos(2*alpha0) , (m*l^2 +m*a^2+m_H*l^2)-m*b*l*cos(2*alpha0) ;
				m*b^2 , -m*b*l*cos(2*alpha0) ];

	disp(['- Is Q(alpha) invertible? Rank[ Q(alpha) ] = ' num2str(rank(Q_alpha)) ] )

	H_alpha = (Q_alpha)^(-1) * P_alpha;

	eps0 = 10^(-2);

	results.Parameters.g = g;
	results.Parameters.M = M; results.Parameters.N = N;	results.Parameters.G = G;

	disp(' ')

	%%%%%%%%%%%%%
	%% Testing %%
	%%%%%%%%%%%%%

	disp('Beginning test 1.')
	phi = pi/6;

	x1 = [ pi+0.5*(pi-2*alpha0)+(pi/2 -phi ) ; alpha0 - phi ; 0.1 ; 0.1 ];
	results.test1.x1 = x1;

	in_transition_region = norm( x1(1:2) - [ pi+0.5*(pi-2*alpha0)+(pi/2 -phi ) ; alpha0 - phi ] ) < eps0;
	results.test1.x1_in_transition_region = in_transition_region;

	x1_prime = x1;
	x1_prime([3:4]) = H_alpha*x1([3:4]);

	results.test1.x1_prime = x1_prime;






end