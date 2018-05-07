classdef FHAE_sb
%Description:
%	This class is the state based controller that was developed in my most recent write up. (#31)
	properties
		L;
		F_set;
		u0_set;
	end

	methods
		%Constructor
		function contr = FHAE_sb( L , F_set , u0_set )
			%Simply copy everything over
			contr.L = L;
			contr.F_set = F_set;
			contr.u0_set = u0_set;
		end

		%Apply the correct control
		function u = apply_control( obj , observed_w , y_vec )
			disp('This function doesn''t work yet!')

			observed_w

			%Useful Constants
			ow_len = length(observed_w);
			m = size(obj.F_set{1},1) / size(obj.L,2);
			p = size(obj.F_set{1},2) / size(obj.L,2);

			%Match observed_w with a word from L
			matching_mat = obj.L(:,1:ow_len);
			matching_mat = repmat(observed_w,size(obj.L,1),1) == matching_mat;

			matching_locs = all(matching_mat')';

			first_match_ind = find(matching_locs,1)

			%Obtain the correct feedback matrices
			obj.F_set{first_match_ind};
			F_t  = obj.F_set{first_match_ind}([m*(ow_len-1)+1:m*ow_len],[1:p*ow_len])
			u0_t = obj.u0_set{first_match_ind}([m*(ow_len-1)+1:m*ow_len])

			u = F_t * y_vec + u0_t;
		end

		%Simulation of a single run
		function x_0_t = simulate_1run( obj , ad , M1 , in_sig )
			%Description:
			%	Uses the information provided in the Affine Dynamics instance ad and
			%	the problem specification M1 along with the given controller
			%	to simulate a single run of the prefix based controller.
			%Usages:
			%
			%Outputs:
			%	x_0_t - An n x (T+1) matrix which defines the trajectory of the state
			%			for a feasible realization of the random variables associated
			%			with the tuple (ad,M1,sigma)

			%Constants
			T = size(obj.L,2);
			n = size(ad.A,1);
			wd = size(ad.B_w,2)
			vd = size(ad.C_v,2);

			%Generate Random Variables
			x0 = unifrnd(-M1,M1,n,1);
			w  = unifrnd(-ad.eta_w,ad.eta_w,wd,T);
			v  = unifrnd(-ad.eta_v,ad.eta_v,vd,T+1);

			if nargin ~= 4
				sig = obj.L( randi(size(obj.L,1),1) , : );
			else
				if (size(in_sig,1) ~= 1) || (size(in_sig,2) ~= T)
					error('Input word is not a single word or does not have the correct length.' )
				end
				sig = in_sig;
			end

			%Simulate system forward.
			x_t = x0;
			y_t = ad.C*x0 + ad.C_v*v(:,1);
			
			x_0_t = x_t;
			y_0_t = y_t;
			x_tp1 = Inf(n,1);
			sig = [sig,1];
			for t = 0:T-1
				%Use Affine Dynamics with proper control law.
				x_tp1 = ad.A * x_t + ...
						ad.B * obj.apply_control( sig([1:t+1]) , y_0_t ) + ...
						ad.B_w * w(:,t+1);
				%Update other variables in system
				x_t = x_tp1;
				y_t = sig(t+2)*(ad.C*x_t + ad.C_v*v(:,t+2));

				x_0_t = [x_0_t x_t];
				y_0_t = [y_0_t; y_t];
			end

		end

		function [x_tensor,x_norm_tensor] = simulate_n_runs( obj , ad , M1 , num_runs , in_sig )
			%Description:
			%	Uses the information provided in the affine dynamics instance ad and the problem
			%	specification M1 along with the given controller to simulate as many runs as the user
			%	would like.
			%
			%Inputs:
			%	x_tensor = simulate_n_runs( obj , ad , M1 , num_runs , in_sig )
			%
			%Outputs:
			%	x_tensor - An n x (T+1) x num_runs tensor which defines the num_runs trajectories of the
			%			   state.

			if nargin < 3
				error('Not enough inputs.')
			end

			%% Constants
			T = size(obj.L,2);

			x_tensor = [];
			for run_ind = 1:num_runs
				%
				if nargin < 5
					x_tensor(:,:,run_ind) = obj.simulate_1run( ad , M1 );
				else
					x_tensor(:,:,run_ind) = obj.simulate_1run( ad , M1 , in_sig );
				end

				%Calculate Norms
				for t = 0 : T
					x_norm_tensor(t+1,:,run_ind) = norm( x_tensor(:,t+1,run_ind) , Inf );
				end
				x_norm_tensor(t+1,:,run_ind)
			end

		end

	end

end
