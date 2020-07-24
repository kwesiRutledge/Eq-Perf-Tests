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
					-m*l*b*cos(theta_s - theta_ns) , (m_H+m)*l^2+m*a^2 ]

		end

		function plot( cw )
			%Description:
			%Plots the compass walker into the current matlab figure.

			%% Constants

			%Ramp Relevant points
			ramp_bounds = [-3,3];

			ramp_top_point = [ramp_bounds(1);0];
			rot_mat = [ cos( -cw.phi ) , -sin( -cw.phi ) ;
						sin( -cw.phi ) , cos( -cw.phi ) ];

			ramp_top_point = rot_mat * [ ramp_bounds(1) - ramp_bounds(2) ; 0 ];
			ramp_top_point = ramp_top_point + [ ramp_bounds(2) ; 0 ];

			%Stationary Leg
			s_leg_point_on_ramp = rot_mat* [ 2 - ramp_bounds(2) ; 0 ];
			s_leg_point_on_ramp = ramp_top_point + [ramp_bounds(2) ; 0];



			%% Algorithm
			hold on;

			%Draw the ramp.
			plot([max(ramp_bounds(1),ramp_top_point(1)),ramp_bounds(2)],zeros(1,2),'b'); %The level ground
			plot([ramp_top_point(1);ramp_bounds(2)],[ramp_top_point(2);0],'b')
			plot(max(ramp_bounds(1),ramp_top_point(1))*ones(2,1),[ramp_top_point(2);0],'b');



		end

	end

end