classdef FHAE_pb
%Description:
%	This class is the prefix-based controller that was developed in my most recent write up. (#31)
	properties
		L;
		F_set;
		u0_set;
	end

	methods
		%Constructor
		function contr = FHAE_pb( L , F_set , u0_set )
			%Description:
			%	Simply copy everything over
			%Inputs:
			%	L - 	A cell array of one dimension.
			%	F_set - A 1 x |L| cell array of feedback matrices which are indexed
			%			based on the index of the matching word in L.
			
			if isnumeric(L)
				%Convert L to cell array.
				L_temp = {};
				for word_ind = 1:size(L,1)
					L_temp{word_ind} = L(word_ind,:);
				end
				contr.L = L_temp;
			elseif iscell(L)
				contr.L = L;
			else
				error('Unknown type of L.')
			end

			contr.F_set = F_set;
			contr.u0_set = u0_set;
		end

		%Apply the correct control
		function u = apply_control( obj , observed_w , y_vec )
			%Useful Constants
			ow_len = length(observed_w);
			num_words = length(obj.L);
			m = size(obj.F_set{1},1) / length(obj.L{1});
			p = size(obj.F_set{1},2) / length(obj.L{1});

			%Match observed_w with a word from L
			matching_mat = [];
			for word_ind = 1:length(obj.L)
				if length(obj.L{word_ind}) >= ow_len
					matching_mat = [matching_mat; obj.L{word_ind}(1:ow_len)]
				else
					%If the word is not long enough to match, then place some nonsense in the
					%corresponding row of the matching matrix.
					matching_mat = [matching_mat; Inf(1,ow_len)]
				end
			end
			matching_mat = repmat(observed_w,num_words,1) == matching_mat;

			matching_locs = all(matching_mat')';

			first_match_ind = find(matching_locs,1);

			%Obtain the correct feedback matrices
			obj.F_set{first_match_ind};
			F_t  = obj.F_set{first_match_ind}([m*(ow_len-1)+1:m*ow_len],[1:p*ow_len]);
			u0_t = obj.u0_set{first_match_ind}([m*(ow_len-1)+1:m*ow_len]);

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
			% T = size(obj.L,2);
			n = size(ad.A,1);
			wd = size(ad.B_w,2);
			vd = size(ad.C_v,2);

			if nargin ~= 4
				rand_word_ind = randi(length(obj.L),1);
				sig = obj.L{rand_word_ind};
				T = length(obj.L{rand_word_ind});
			else
				if (size(in_sig,1) ~= 1) %|| (size(in_sig,2) ~= T)
					error('Input word is not a single word or does not have the correct length.' )
				end
				sig = in_sig;
			end

			%Generate Random Variables
			x0 = unifrnd(-M1,M1,n,1);
			w  = unifrnd(-ad.eta_w,ad.eta_w,wd,T);
			v  = unifrnd(-ad.eta_v,ad.eta_v,vd,T);

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
						ad.B_w * w(:,t+1) + ...
						ad.f ;
				%Update other variables in system
				x_t = x_tp1;
				x_0_t = [x_0_t x_t];

				if( t == T-1 )
					continue;
				end

				y_t = sig(t+2)*(ad.C*x_t + ad.C_v*v(:,t+1));
				y_0_t = [y_0_t; y_t];
			end

		end

		function [run_data_x,run_data_x_norm] = simulate_n_runs( obj , ad , M1 , num_runs , in_sig )
			%Description:
			%	Uses the information provided in the affine dynamics instance ad and the problem
			%	specification M1 along with the given controller to simulate as many runs as the user
			%	would like.
			%
			%Inputs:
			%	x_tensor = simulate_n_runs( obj , ad , M1 , num_runs , in_sig )
			%
			%Outputs:
			%	run_data_x -	A 1 x 'num_runs' cell array of n x (T_i+1) matrices which defines
			%					the num_runs trajectories of the state.
			%			   		T_i can change for each run.
			%	run_data_x_norm - 	A cell array of

			if nargin < 3
				error('Not enough inputs.')
			end

			%% Constants
			n = size(ad.A,1);

			%% Implement Loop
			run_data_x = {};
			run_data_x_norm = {};
			for run_ind = 1:num_runs
				%
				if nargin < 5
					run_data_x{run_ind} = obj.simulate_1run( ad , M1 );
				else
					run_data_x{run_ind} = obj.simulate_1run( ad , M1 , in_sig );
				end

				%Calculate Norms
				T = size(run_data_x{run_ind},2)-1;
				for t = 0 : T
					run_data_x_norm{run_ind}(t+1,:) = norm( run_data_x{run_ind}(:,t+1) , Inf );
				end
			end

		end

	end

end
