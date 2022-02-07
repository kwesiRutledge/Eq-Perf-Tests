classdef InvertedThickPendulum
	%Description:
	%	Inverted Thick Pendulum

	properties
		CoMx_rel;
		CoMy_rel;
		Height;
		Length;
		Mass;
		mu_rot; %Rotational Coefficient of Friction

		%State
		x; 	% Two dimensional vector. First element is the angle 
			% between the vertical and the ray pointing 
			% towards the upper left hand side of the pendulum. (according to drawing below)
	end

	methods
		function systemOut = InvertedThickPendulum(varargin)
			%InvertedThickPendulum
			%Description:
			%	Initializer for the Inverted Thick Pendulum
			%
			%	  |<---- L --->|
			%	_ ______________ _
			%	^ |            | ^
			%	| |            | |
			%	| |            | Height
			%	| |            | |
			%	v |            | V
			%	_ ______________ _
			%
			%	Where the axis of rotation is on the bottom right corner.
			%
			%	The center of mass is located at the point (CoMx_rel,CoMy_rel) relative
			%	to the reference frame starting in the BOTTOM RIGHT corner of the square.
			%
			%Usage:
			%	itp = InvertedThickPendulum()
			%

			%% Input Processing

			% varargin_idx = 1;
			% while varargin_idx <= nargin
			% 	switch varargin{varargin_idx}
			% 		case 'm'
			% 			body
			% 		otherwise
			% 			body
			% 	end
			% end

			%% Defaults

			if ~isfield('Mass',systemOut)
				systemOut.Mass = 1;
			end

			systemOut.Length = 1.5;
			systemOut.Height = 1.0;

			systemOut.CoMx_rel = (systemOut.Length/2) - systemOut.Length; %Relative to the pendulum's frame.
			systemOut.CoMy_rel = (systemOut.Height/2);  

			systemOut.mu_rot = 0; %Set coefficient of friction to be zero by default.

			if ~isfield(systemOut,'x')
				systemOut.x = [0.1;0.05];
			end

		end

		function [plot_handles] = Show( itp )
			%Description:
			%	Plots the current state of the Inverted Thick Pendulum.
			%	Highlights the points that 

			%% Constants

			lr_corner_pos = [0;0];

			alpha_ang = atan(itp.Height/itp.Length);
			theta = wrapToPi(itp.x(1));

			lw = 2.0; %LineWidth for Rectangle/Pendulum
            pendulumColorChoice = 'blue';

            markerSizeCoM = 72; %Size of the Marker for Center Of Mass Position
            CoMColorChoice = 'blue';

			%% Construct All Geometric Features

			% Plotting All Corners of the Thick Pendulum
			corners = [ lr_corner_pos(1) , lr_corner_pos(2) ;
						lr_corner_pos(1) , lr_corner_pos(2) + itp.Height ;
						lr_corner_pos(1) - itp.Length , lr_corner_pos(2) + itp.Height ;
                        lr_corner_pos(1) - itp.Length , lr_corner_pos(2) ]';

            if (0 <= theta) && (theta <= pi)
            	angle_wrt_horizontal = -((pi/2) - theta - alpha_ang) ;
            elseif (-pi <= theta)  && (theta <= 0)
            	angle_wrt_horizontal = -(pi/2 + (-theta ) - alpha_ang );
            else
            	error(['Unexpected angle value ' num2str(theta)])
            end
            	

            rot_phi = [ cos( angle_wrt_horizontal ) , -sin( angle_wrt_horizontal ) ;
            			sin( angle_wrt_horizontal ) , cos( angle_wrt_horizontal ) ];

            rot_corners = rot_phi * corners;

            %% Plot

            hold on;
            plot_handles = [];

            % Plot Sides
            for corner_idx = 1:size(rot_corners,2)-1
                plot_handles(corner_idx) = plot( ...
                    rot_corners(1,corner_idx+[0:1]) , ...
                    rot_corners(2,corner_idx+[0:1]) , ...
                    'LineWidth', lw , ...
                    'Color', pendulumColorChoice ...
                );
            end
            plot_handles(length(plot_handles)+1) = plot( ...
                    rot_corners(1,[4,1]) , ...
                    rot_corners(2,[4,1]) , ...
                    'LineWidth', lw , ...
                    'Color', pendulumColorChoice ...
            );

            % Plot Center of Mass
            CoM_coords = [ itp.CoMx_rel ; itp.CoMy_rel ];
            rot_CoM = rot_phi * CoM_coords;

            scatter( ...
            	rot_CoM(1),rot_CoM(2), ...
            	markerSizeCoM , CoMColorChoice , 's' , 'filled' ...
            )

            % Plot the Upper Left Corner (represented by state)
            scatter( ...
            	rot_corners(1,3), rot_corners(2,3), ...
            	markerSizeCoM , CoMColorChoice , 'd' , 'filled' ...
            )

            hold off;

		end

		function theta_tilde = get_CoM_angle_offset(itp)
			%Description:
			%	The position of the center of mass causes there to be an offset between
			%	theta (x(1)) and the true angle at which the center of mass is located
			%	theta + theta_tilde. This finds the offset theta_tilde.

			%% Constants

			%% Algorithm

			theta_CoM0 = atan( itp.CoMy_rel / abs(itp.CoMx_rel) );

			theta_tilde = - (theta_CoM0 - atan( itp.Height / itp.Length ));	

		end

		function dxdt = f( itp , x , u )
			%Description:
			%	Describes what the derivative of the system's dynamics are for a given
			%	state and input u.

			%% Constants

			m = itp.Mass;
			g = 9.8;

			r_CoM = norm( [itp.CoMx_rel ; itp.CoMy_rel] ,2);

			%% Algorithm

			dxdt(1) = x(2);
			dxdt(2) = r_CoM * m * g * cos( (pi/2) - (x(1) + itp.get_CoM_angle_offset()) ) - itp.mu_rot * x(2) + u;

			dxdt = dxdt'; %Make this a column vector.

		end

		function [Ac,Bc,Kc] = LinearizedContinuousDynamicsAbout( itp , x , u )
			%Description:
			%	Creates the linearized, continuous time system when the linearization is centered at the state x and input u
			%	and evaluated at the current state.
			%
			%Usage:
			%	[ Ac , Bc ] = itp.LinearizedContinuousDynamicsAbout( x , u )

			% Constants
			syms theta_x theta_dot_x symvar_u real;

            symvar_x = [theta_x; theta_dot_x];

            % Dynamics
            dxdt = itp.f(symvar_x,symvar_u);

            dfdx = jacobian(dxdt,symvar_x);
            dfdu = jacobian(dxdt,symvar_u);

            %Give the symbolic variables real values based
            %on the input x and u
            theta_x = x(1); theta_dot_x = x(2);
            symvar_u = u(1); 

            Ac = eval( subs(dfdx) );
            Bc = eval( subs(dfdu) );

           	Kc = itp.f( x , u );

        end

        function [ Ad , Bd , Kd ] = LinearizedDiscreteDynamicsAbout( itp , x , u , dt )
			%Description:
			%	Creates the linearized, discrete-time system when the linearization is centered at the state x and input u
			%	and uses discretization state.

			% Constants

            % Get Continuous Time Matrices
            [ Ac , Bc , Kc ] = itp.LinearizedContinuousDynamicsAbout( x , u );

            n = size(Ac,1);

            % Discretize using MATLAB's built-in tools
            csys1 = ss(Ac,Bc,eye(n),0);
            dsys1 = c2d(csys1,dt);

            csys2 = ss(Ac,Kc,eye(n),0);
            dsys2 = c2d(csys2,dt);

            Ad = dsys1.A;
            Bd = dsys1.B;
            Kd = dsys2.B;

        end

	end

end