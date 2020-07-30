classdef CompassWalker
	%CompassWalker Class
	%Description:
	%

	properties(Constant)
		g = 10;
	end

	properties
		a; 		%Distance along leg to the Center of Mass
		b; 		%Ditance from leg Center of Mass to Center of Mass of Walker's "body"
		m; 		%Mass of each leg
		m_H; 	%Mass of the walker's "body"
		l; 		%Length of the leg
		phi; 	%Angle of Walking Ramp
		CurrentState; 	%Current state
	end

	methods
		function cw = CompassWalker( varargin )
			%CompassWalker.m
			%Description:
			%	This constructor allows one to define a Compass Walker object with:
			%	- 
			%
			%Usage:
			%	cw = CompassWalker()
			%	cw = CompassWalker('WalkerConstants',a,b,m,m_H)
			%	cw = CompassWalker('CurrentState')

			%% Input Processing %%
			
			argidx = 1;
			while nargin >= argidx
				switch varargin{argidx}
					case 'WalkerConstants'
						cw.a = varargin{argidx+1};
						cw.b = varargin{argidx+2};
						cw.m = varargin{argidx+3};
						cw.m_H = varargin{argidx+4};

						argidx = argidx + 5;
					case 'RampAngle'
						cw.phi = varargin{argidx+1};
						argidx = argidx + 2;
					case 'CurrentState'
						cw.CurrentState = varargin{argidx+1};
						argidx = argidx + 2;
					otherwise
						error('Unrecognized input to CompassWalker().')
				end
			end

			%% Constants %%

			if ~isfield(cw,'a')
				cw.a = 0.5; 	%meters
			end

			if ~isfield(cw,'b')
				cw.b = 0.5;
			end

			if ~isfield(cw,'m')
				cw.m = 5;
			end

			if ~isfield(cw,'m_H')
				cw.m_H = 10;
			end

			cw.l = cw.a+cw.b;

			if ~isfield(cw,'phi')
				cw.phi = pi/6; %30 degrees, or pi/6 radians
			end

			if ~isfield(cw,'CurrentState')
				cw.CurrentState = NaN(4,1);
			end

			%% Algorithm %%

		end

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%% Defining State-Dependent Dynamics %%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		function M_out = M( cw , theta_in )
			%Description:
			%
			%Inputs:
			%	theta_in - 2 x 1 input vector representing 
			%				[ theta_ns ]
			%				[ theta_s  ]

			%% Constants %%

			m = cw.m;
			m_H = cw.m_H;
			a = cw.a;
			b = cw.b;
			l = cw.l;

			%% Algorithm %%

			theta_ns = theta_in(1);
			theta_s = theta_in(2);

			M_out = [ m*(b^2) , -m*l*b*cos( theta_s - theta_ns ) ;
					-m*l*b*cos(theta_s - theta_ns) , (m_H+m)*l^2+m*a^2 ];

		end

		function N_out = N( cw , x )
			%Description:
			%
			%Usage:
			%	N_out = cw.N()
			%	N_out = cw.N( theta )
			%

			% %% Input Checking

			% cw = varargin{1};

			% argidx = 2;
			% if nargin == 2
			% 	x = varargin{2};
			% end

			%% Constants

			% if ~exist('x')
			% 	if all(isnan(cw.CurrentState))
			% 		error('The Current State of the Compass Walker is not defined!')
			% 	elseif isnumeric(cw.CurrentState)
			% 		x = cw.CurrentState;
			% 	else
			% 		error(['CurrentState is an unexpected type!']);
			% 	end
			% end

			theta_ns 	 = x(1);
			theta_s 	 = x(2);
			dot_theta_ns = x(3);
			dot_theta_s  = x(4);

			m = cw.m;
			m_H = cw.m_H;
			a = cw.a;
			b = cw.b;
			l = cw.l;

			%% Algorithm

			N_out = [ 	0 , m*l*b*dot_theta_s*sin(theta_s - theta_ns);
					-m*l*b*dot_theta_ns*sin(theta_s - theta_ns) , 0 ];

		end

		function G_out = G( cw , x )
			%

			% %% Input Checking

			% cw = varargin{1};

			% argidx = 2;
			% if nargin == 2
			% 	x = varargin{2};
			% end

			%% Constants

			% if ~exist('x')
			% 	if all(isnan(cw.CurrentState))
			% 		error('The Current State of the Compass Walker is not defined!')
			% 	elseif isnumeric(cw.CurrentState)
			% 		x = cw.CurrentState;
			% 	else
			% 		error(['CurrentState is an unexpected type!']);
			% 	end
			% end

			theta_ns 	 = x(1);
			theta_s 	 = x(2);
			dot_theta_ns = x(3);
			dot_theta_s  = x(4);

			m = cw.m;
			m_H = cw.m_H;
			a = cw.a;
			b = cw.b;
			l = cw.l;
			g = cw.g;

			%% Algorithm

			G_out = [ m*b*sin(theta_ns) ; -(m_H*l+m*a+m*l)*sin(theta_s) ] * g;

		end
			

		function dxdt = cDynamics(varargin)
			%Description:
			%	Defines the derivative of the Compass Walker's continuous dynamics
			%	when the current time and state (x) are given.
			%
			%Usage:
			%	dxdt = cw.cDynamics(u)
			%	dxdt = cw.cDynamics(u,'x',x)

			%% Input Processing
			u = varargin{1};

			if any(size(u) ~= [3,1])
				error(['Input must be a 3x1 vector. Instead, it is ' num2str(size(u,1)) 'x' num2str(size(u,2)) '.' ])
			end

			argidx = 2;
			while argidx <= nargin
				switch varargin{argidx}
					case 'x'
						x = varargin{argidx+1};
						argidx = argidx + 2;
					otherwise
						error('Unexpected input to cDynamics().')
				end
			end

			%% Constants

			if ~exist('CurrentState')
				if all(isnan(cw.CurrentState))
					error('The Current State of the Compass Walker is not defined!')
				elseif isnumeric(cw.CurrentState)
					x = cw.CurrentState;
				else
					error(['CurrentState is an unexpected type!']);
				end
			end

			S = [	1,-1,0;
					0,1,1];

			%% Algorithm

			dxdt = [	x(3) ; x(4) ;
						cw.M([x(1);x(2)])^(-1) * ...
						( -cw.N([x(1),x(2),x(3),x(4)]')*[x(3);x(4)] + S*u - cw.G(x(1),x(2)) ) ];

		end


		%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%% Visualizing Functions %%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%

		function [plot_handles] = plot( cw )
			%Description:
			%Plots the compass walker into the current matlab figure.

			%% Constants

			rot_mat = @(angle_in) [ cos( angle_in ) , -sin( angle_in ) ;
									sin( angle_in ) , cos( angle_in ) ];

			%Ramp Relevant points
			ramp_bounds = [-3,3];

			ramp_top_point = [ramp_bounds(1);0];
			rot_mat_phi = rot_mat(-cw.phi);

			ramp_top_point = rot_mat_phi * [ ramp_bounds(1) - ramp_bounds(2) ; 0 ];
			ramp_top_point = ramp_top_point + [ ramp_bounds(2) ; 0 ];

			ramp_bottom_point = [3;0];

			%Stationary Leg
			alpha_s = 0.65;
			s_leg_point_on_ramp = alpha_s*ramp_bottom_point + (1-alpha_s)*ramp_top_point;			

			s_leg_pt1 = [ 0 ; 0 ];
			s_leg_pt2 = [ cw.l ; 0 ];

			s_leg_pt2 = rot_mat(pi/2) * rot_mat(cw.CurrentState(2)) * s_leg_pt2;

			s_leg_pt1 = s_leg_point_on_ramp + s_leg_pt1;
			s_leg_pt2 = s_leg_point_on_ramp + s_leg_pt2;
			
			body_pos = s_leg_pt2;

			%Plot Body Position
			body_drawing.r = 0.1;
			body_drawing.center = body_pos;
			body_drawing.pos = [ 	body_drawing.center' - body_drawing.r , ...
									2*body_drawing.r , 2*body_drawing.r ];

			%Nonstationary Leg
			ns_leg_pt1 = body_pos;
			ns_leg_pt2 = [ -cw.l ; 0 ];

			ns_leg_pt2 = rot_mat( pi/2 - (2*pi-cw.CurrentState(1)) ) * ns_leg_pt2;
			ns_leg_pt2 = body_pos + ns_leg_pt2;


			%% Algorithm

			%Draw the ramp.
			plot_handles(1) = plot([max(ramp_bounds(1),ramp_top_point(1)),ramp_bounds(2)],zeros(1,2),'b'); %The level ground
			hold on;
			plot_handles(2) = plot([ramp_top_point(1);ramp_bounds(2)],[ramp_top_point(2);0],'b');
			plot_handles(3) = plot(max(ramp_bounds(1),ramp_top_point(1))*ones(2,1),[ramp_top_point(2);0],'b');

			%Plot Stationary leg
			plot_handles(4) = plot([ s_leg_pt1(1) , s_leg_pt2(1) ],[ s_leg_pt1(2) , s_leg_pt2(2) ],'r');

			%Plot Body Position
			plot_handles(5) = rectangle('Position',body_drawing.pos,'Curvature',[1,1]);

			%Plot Stationary leg
			plot_handles(6) = plot([ ns_leg_pt1(1) , ns_leg_pt2(1) ],[ ns_leg_pt1(2) , ns_leg_pt2(2) ],'g');


			axis equal

		end

		function visualize_trajectory( varargin )
			%Description:
			%	Simulates the given trajectory of states for the compass walker.
			%
			%Usage:
			%	cw.visualize_trajectory( x_trajectory )
			%	cw.visualize_trajectory( x_trajectory , 'FileName' , '../my_own_video.mp4' )
			%	cw.visualize_trajectory( x_trajectory , 'StoppingIndex' , 200 )

			%% Input Checking
			cw = varargin{1};
			x_trajectory = varargin{2};

			if size(x_trajectory,1) ~= 4
				%
				error(['The input trajectory x_trajectory should have 4 rows, instead it has ' num2str(size(x_trajectory,1)) '.' ])
			end

			argidx = 3;
			while argidx <= nargin
				switch varargin{argidx}
					case 'FileName'
						FileName = varargin{argidx+1};
						argidx = argidx + 2;
					case 'StoppingIndex'
						StoppingIndex = varargin{argidx + 1};
						argidx = argidx + 2;
					otherwise
						error(['Unexpected input to visualize_trajectory(): ' varargin{argidx} ])
				end
			end

			%% Constants

			traj_length = size(x_trajectory,2);

			if ~exist('FileName')
				FileName = 'results/walker_vid.mp4';
			end

			vidObj = VideoWriter(FileName,'MPEG-4');

			if ~exist('StoppingIndex')
				StoppingIndex = traj_length;
			end

			%% Algorithm

			open(vidObj);

			figure;
			for t = 1:StoppingIndex
				%Collect State from the input list.
				cw.CurrentState = x_trajectory(:,t);

				%Plot
				h = cw.plot();
				
				%Get Current frame and write it.
				currFrame = getframe;
				writeVideo(vidObj,currFrame);
				delete(h);

				%Prepare for plot to be overwritten
				hold off; 

			end

			close(vidObj);

		end

		function alpha_out = alpha(cw)
			%Description:
			%	Computes the interleg angle 2*alpha using the CurrentState
			%	of the walker.

			%% Input Processing

			if isnan(cw.CurrentState)
				error('CurrentState is not initialized yet!')
			end

			%% Constants

			rot_mat = @(angle_in) [ cos( angle_in ) , -sin( angle_in ) ;
									sin( angle_in ) , cos( angle_in ) ];

			%Ramp Relevant points
			ramp_bounds = [-3,3];

			ramp_top_point = [ramp_bounds(1);0];
			rot_mat_phi = rot_mat(-cw.phi);

			ramp_top_point = rot_mat_phi * [ ramp_bounds(1) - ramp_bounds(2) ; 0 ];
			ramp_top_point = ramp_top_point + [ ramp_bounds(2) ; 0 ];

			ramp_bottom_point = [3;0];

			%Stationary Leg
			alpha_s = 0.65;
			s_leg_point_on_ramp = alpha_s*ramp_bottom_point + (1-alpha_s)*ramp_top_point;			

			s_leg_pt1 = [ 0 ; 0 ];
			s_leg_pt2 = [ cw.l ; 0 ];

			s_leg_pt2 = rot_mat(pi/2) * rot_mat(cw.CurrentState(2)) * s_leg_pt2;

			s_leg_pt1 = s_leg_point_on_ramp + s_leg_pt1;
			s_leg_pt2 = s_leg_point_on_ramp + s_leg_pt2;
			
			body_pos = s_leg_pt2;

			%Plot Body Position
			body_drawing.r = 0.1;
			body_drawing.center = body_pos;
			body_drawing.pos = [ 	body_drawing.center' - body_drawing.r , ...
									2*body_drawing.r , 2*body_drawing.r ];

			%Nonstationary Leg
			ns_leg_pt1 = body_pos;
			ns_leg_pt2 = [ -cw.l ; 0 ];

			ns_leg_pt2 = rot_mat( pi/2 - (2*pi-cw.CurrentState(1)) ) * ns_leg_pt2;
			ns_leg_pt2 = body_pos + ns_leg_pt2;

			%% Algorithm

			triangle_side = norm( ns_leg_pt2 - s_leg_point_on_ramp ,2);

			alpha_out = (pi/2) - acos(triangle_side/2*cw.l);

		end


	end

end