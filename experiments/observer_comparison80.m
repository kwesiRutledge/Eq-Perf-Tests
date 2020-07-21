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

	Q_alpha = [ m*b^2 - m*b*l*cos(2*alpha0) , (m*l^2 + m*a^2+m_H*l^2)-m*b*l*cos(2*alpha0) ;
				m*b^2 , -m*b*l*cos(2*alpha0) ];

	P = @(alpha_in) ...
		[ (m_H*l^2 + 2*m*l^2)*cos(2*alpha_in) - m*a*b - 2*m*b*l*cos(2*alpha_in) , -m*a*b ;
			-m*a*b , 0 ];
	Q = @(alpha_in) ...
		[ m*b^2 - m*b*l*cos(2*alpha_in) , (m*l^2 + m*a^2+m_H*l^2)-m*b*l*cos(2*alpha_in) ;
				m*b^2 , -m*b*l*cos(2*alpha_in) ];

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
	disp('Attempting to verify some of the transfer conditions.')

	slope_grade = pi/6; %Current slope. phi
	switching_angle = alpha0; %Distance of current. alpha
	
	results.test1.slope_grade = slope_grade;
	results.test1.switching_angle = switching_angle;

	H_t = ( Q(switching_angle) )^(-1) * P(switching_angle);

	x1 = zeros(4,1);
	%theta1 = (7*pi/4) - switching_angle;

	theta_s = 7*pi/4;
	theta1 = theta_s - ( pi + (pi/2) - slope_grade );

	% if slope_grade >= switching_angle
	% 	theta_s = pi + (pi-theta1) + ( pi/2 - slope_grade );
	% else
	% 	theta_s = pi - (pi/2 - slope_grade ) - theta1;
	% end

	if slope_grade >= switching_angle
		theta_ns = (pi/2) - slope_grade - theta1;
	else
		theta_ns = pi + (pi - theta1) + (pi/2) - slope_grade;
	end

	

	x1(1) = theta_ns; %theta_ns
	x1(2) = theta_s; %theta_s
	x1(3) = 0.1; %theta_ns - dot
	x1(4) = -0.1; %theta_s - dot

	results.test1.x1 = x1;

	abs( wrapToPi(x1(1))+wrapToPi(x1(2))+2*slope_grade)

	tf_ns_s_relationship = abs( wrapToPi(x1(1)) + wrapToPi(x1(2)) + 2*slope_grade) < eps0;
	disp(['- x1(1) + x1(2) == -2slope_grade? ' num2str(tf_ns_s_relationship) ])

	in_transition_region = norm( x1(1:2) - [ pi+0.5*(pi-2*alpha0)+(pi/2 -slope_grade ) ; alpha0 - slope_grade ] ) < eps0;
	disp(['- in_transition_region = ' num2str(in_transition_region) ])

	x1_prime = x1;
	x1_prime([3:4]) = H_t*x1([3:4]);

	disp('- x_t = ')
	x1
	disp('- x1_prime = ')
	x1_prime

	results.test1.x1_prime = x1_prime;
	results.test1.x1_in_transition_region = in_transition_region;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Test 2: Swing Leg Trajectory %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	disp('Test 2: Swing Leg Trajectory')
	disp('')

	x2 = zeros(4,1);
	x2(1) = 7*pi/4; %Initial theta_ns
	x2(2) = 11*pi/6; %Initial theta_s
	x2(3) = 0;
	x2(4) = 0;

	u2_constant = zeros(3,1);

	dynamics = @(t,x) ...
				[	x(3) ; x(4) ;
					M(x(1),x(2))^(-1) * ...
					( -N(x(1),x(2),x(3),x(4))*[x(3);x(4)] + S*u2_constant - G(x(1),x(2)) ) ];

	tspan2 = [0,1];

	[t2_out,x2_out] = ode45(dynamics,tspan2,x2);

	%Find the point at which 
	ns_leg_passes_s_leg = x2_out(:,1) < x2_out(:,2);

	disp(['- Total Number of Points in Trajectory = ' num2str(length(t2_out)) ])
	disp(['- Number of points where non-stationary theta is larger than stationary theta = ' num2str(sum(ns_leg_passes_s_leg)) ])

	%Plotting
	figure;
	for state_idx = 1:size(x2_out,2)
		subplot(4,1,state_idx)
		plot(t2_out,x2_out(:,state_idx))
	end

	disp(' ')

	results.test2.x2 = x2;
	results.test2.x2_out = x2_out;
	results.test2.ns_leg_passes_s_leg = ns_leg_passes_s_leg;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Test 3: Swing Leg Trajectory 2 %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	disp('Test 3: Swing Leg Trajectory with inputs')
	disp('')

	x3 = zeros(4,1);
	x3(1) = 7*pi/4; %Initial theta_ns
	x3(2) = 11*pi/6; %Initial theta_s
	x3(3) = 0;
	x3(4) = 0;

	u3_constant = zeros(3,1);

	dynamics3 = @(t,x,u) ...
				[	x(3) ; x(4) ;
					M(x(1),x(2))^(-1) * ...
					( -N(x(1),x(2),x(3),x(4))*[x(3);x(4)] + S*u - G(x(1),x(2)) ) ];

	tspan3 = [0,1];

	[t3_out,x3_out] = ode45(@(t,x) dynamics3(t,x,u3_constant),tspan3,x3);

	%Find the point at which 
	ns_leg_passes_s_leg = x3_out(:,1) < x3_out(:,2);

	disp(['- Total Number of Points in Trajectory = ' num2str(length(t3_out)) ])
	disp(['- Number of points where non-stationary theta is larger than stationary theta = ' num2str(sum(ns_leg_passes_s_leg)) ])

	%Plotting
	figure;
	for state_idx = 1:size(x3_out,2)
		subplot(4,1,state_idx)
		plot(t3_out,x3_out(:,state_idx))
	end

	disp(' ')

	results.test3.x3 = x3;
	results.test3.x3_out = x3_out;
	results.test3.ns_leg_passes_s_leg = ns_leg_passes_s_leg;

	%%%%%%%%%%%%%%%%%%%
	%% Test 4: Class %%
	%%%%%%%%%%%%%%%%%%%

	disp('Test 4: Creating a Class which will handle a lot of the grunt work.')

	cw4 = CompassWalker();

	x4 = zeros(4,1);
	x4(1) = 7*pi/4; %Initial theta_ns
	x4(2) = 11*pi/6; %Initial theta_s
	x4(3) = 0;
	x4(4) = 0; 

	assert( all(all( M(x4(1),x4(2)) == cw4.M(x4(1:2)) )) ) 

end