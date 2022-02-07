%Quadrotor0.m
classdef Quadrotor0
	%Description:
	%	Quadrotor's dynamics.
	%	With values of constants that match the CrazyFlie drones in the lab.

	properties
		Mass;
		MomentOfIntertia;

		Ixx;
		Iyy;
		Izz;

		% State
		x; 	% Twelve dimensional vector.
	end

	methods
		function systemOut = Quadrotor0(varargin)
			%Quadrotor0
			%Description:
			%
			%
			%Usage:
			%	itp = InvertedThickPendulum()
			%

			%% Input Processing



			%% Defaults

			systemOut.Mass = 0.03; %kg
			
			systemOut.Ixx = 1.395e-5;
			systemOut.Iyy = 1.436e-5;
			systemOut.Izz = 2.173e-5;

			systemOut.MomentOfIntertia = diag([ systemOut.Ixx ; systemOut.Iyy ; systemOut.Izz ]);

			if ~isfield(systemOut,'x')
				systemOut.x = zeros(12,1);
			end

		end

		function [ x_dot ] = f(qr,x,u)
			%Description:
			%	Nonlinear dynamics of the quadrotor.

			% Constants

			r_x = x(1); r_y = x(2); r_z = x(3);
			alpha = x(4); beta = x(5); gamma = x(6);
			v_x = x(7); v_y = x(8); v_z = x(9);
			alpha_dot = x(10); beta_dot = x(11); gamma_dot = x(12);

			m = qr.Mass;

			g = 9.81; % m/s^2

			I_x = qr.Ixx;
			I_y = qr.Iyy;
			I_z = qr.Izz;

			% Algorithm

			x_dot = [
				v_x ;
				v_y ;
				v_z ;
				beta_dot * (sin(gamma)/cos(beta)) + gamma_dot * (cos(gamma)/cos(beta)) ;
				beta_dot * cos(gamma) - gamma_dot * sin(gamma) ;
	        	alpha_dot + beta_dot * sin(gamma) * tan(beta) + gamma_dot * cos(gamma) * tan(beta) ;
	        	-(1/m) * ( sin(gamma) * sin(alpha) + cos(gamma) * cos(alpha) * sin(beta) ) * u(1) ;
		        -(1/m) * ( sin(gamma) * cos(alpha) - cos(gamma) * sin(alpha) * sin(beta) ) * u(1) ;
		        g - (1/m) * cos(gamma) * cos(beta) * u(1) ;
		        ((I_y - I_z)/I_x) * beta_dot * gamma_dot + (1/I_x) * u(2) ;
		        ((I_z - I_x)/I_y) * alpha_dot * gamma_dot + (1/I_y) * u(3) ;
		        ((I_x - I_y)/I_z) * alpha_dot * beta_dot + (1/I_z) * u(4)
			];

		end

		function [ ad0 ] = ToAffineDynamicsLinearizedAbout(qr,x,u)
			%Description:
			%	Linearizes the dynamics about the state x with input u.
			%	Also adds the constant terms in the dynamics which should make the model affine.



		end

	end

end