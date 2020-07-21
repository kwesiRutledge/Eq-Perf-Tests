classdef CompassWalker
	%CompassWalker Class
	%Description:
	%

	properties(Constant)
		g = 10;
	end

	properties
		a;
		b;
		m;
		m_H;
		l;
	end

	methods
		function cw = CompassWalker( varargin )
			%CompassWalker.m
			%Description:
			%

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

	end

end