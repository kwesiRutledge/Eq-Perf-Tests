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
			obj.F_set{first_match_ind}
			F_t  = obj.F_set{first_match_ind}([m*(ow_len-1)+1:m*ow_len],[1:p*ow_len])
			u0_t = obj.u0_set{first_match_ind}([m*(ow_len-1)+1:m*ow_len])

			u = F_t * y_vec + u0_t;
		end

		%Simulation of a single run
		function x_t = simulate_1run( obj , ad , M1 )
			%Description:
			%	Uses the information provided in the Affine Dynamics instance ad and
			%	the problem specification M1 along with the given controller
			%	to simulate a single run of the state based controller.

			%Constants
			T = size(obj.L,2);
			n = size(ad.A,1);
			wd = size(ad.B_w,2);
			vd = size(ad.C_v,2);

			%Generate Random Variables
			x0 = rand(-M1,M1,n,1);
			w  = rand(-ad.eta_w,ad.eta_w,T*wd,1);
			v  = rand(-ad.eta_v,ad.eta_v,T*vd,1);

			%Simulate system forward.
			x_t = x0;
			x_tp1 = Inf(n,1);
			for t = 0:T-1


		end

	end

end
