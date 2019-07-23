function [results] = observer_comparison47(varargin)
	%observer_comparison47.m
	%
	%Description:
	%	The objective of this is to prototype a form of "leader-follower" dynamics. Version 4.
	%
	%Assumptions:
	%	- I am currently assuming that the value of the initial estimation error and the state are the same.
	%	  This only works if we assume that the initial estimate is the zero vector.

	%% Constants
	%Model Parameters
	h = 0.2;
	num_sys = 5;
	tau_i = [0.1, 0.5,0.5,2,2];

	n_i = 3;
	p_i = n_i;

	eta_w_i = [0.06,0.0015,0.0015,0.0005,0.0005];
	eta_v = 0.002;

	%Meta Model Parameters
	r = 2; %Dimension for the follower cube.

	% ad_arr = [];
	for sys_idx = 1:num_sys
		temp_A = [1,h,0;0,1,h; 0,0,1-h/tau_i(sys_idx)];
		temp_B = [0;0;h/tau_i(sys_idx)];
		temp_B_w = [0;0;1];

		ad_arr(sys_idx) = Aff_Dyn(temp_A,temp_B,zeros(n_i,1),eye(n_i), eta_w_i(sys_idx), eta_v, temp_B_w , eye(p_i) );
		% ad_arr(sys_idx) = Aff_Dyn(kron(eye(2),temp_A),kron(eye(2),temp_B),zeros(2*n_i,1),eye(2*n_i), eta_w_i(sys_idx), eta_v, kron(eye(2),temp_B_w) , eye(2*p_i) );

	end

	results.ad_arr = ad_arr;

	%% Generate Reference Trajectories

	contr_durs = [4,4,4];
	contr_vals = [1,0.5,0.75];

	x_r = zeros(n,1,num_sys);

	for time_idx = 1:sum(contr_durs)
		for dur_idx = 1:length(contr_durs)
			if dur_idx == 1
				if (time_idx >= 1) && (time_idx <= sum(contr_durs)
					%
					x_k = x_r(:,time_idx+1,sys_idx);
					x_kp1 = ad_arr(sys_idx).A*x_k + ad_arr(sys_idx).B*contr_vals(dur_idx);
					x_r(:,time_idx+1,sys_idx) = x_kp1;
				end
			else
				if (time_idx > sum(contr_durs(1:dur_idx-1))) && (time_idx <= sum(contr_durs(1:dur_idx)))
					x_k = x_r(:,time_idx+1,sys_idx);
					x_kp1 = ad_arr(sys_idx).A*x_k + ad_arr(sys_idx).B*contr_vals(dur_idx);
					x_r(:,time_idx+1,sys_idx) = x_kp1;
				end
			end
		end
	end


end