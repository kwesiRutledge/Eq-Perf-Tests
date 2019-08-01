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
					matching_mat = [matching_mat; obj.L{word_ind}(1:ow_len)];
				else
					%If the word is not long enough to match, then place some nonsense in the
					%corresponding row of the matching matrix.
					matching_mat = [matching_mat; Inf(1,ow_len)];
				end
			end
			matching_mat = repmat(observed_w,num_words,1) == matching_mat;	
            
            %Showing which words match the observed prefix
            if ow_len == 1
                matching_locs = matching_mat;
            elseif ow_len > 1
                matching_locs = all(matching_mat')';
            else
                error('ow_len is a nonpositive integer.')
            end

			first_match_ind = find(matching_locs,1);

			%Obtain the correct feedback matrices
			obj.F_set{first_match_ind};
			F_t  = obj.F_set{first_match_ind}([m*(ow_len-1)+1:m*ow_len],[1:p*ow_len]);
			u0_t = obj.u0_set{first_match_ind}([m*(ow_len-1)+1:m*ow_len]);

			u = F_t * y_vec + u0_t;
		end

		%Simulation of a single run
		function x_0_t = simulate_1run( varargin )
			%Description:
			%	Uses the information provided in the Affine Dynamics instance ad and
			%	the problem specification M1 along with the given controller
			%	to simulate a single run of the prefix based controller.
			%Usages:
			%	x_0_t = obj.simulate_1run( ad , M1 , in_sig )
			%	x_0_t = obj.simulate_1run( ad , M1 , 'in_sigma' , in_sig )
			%	x_0_t = obj.simulate_1run( ad , M1 , 'in_sigma' , in_sig , 'in_w' , in_w )
			%	x_0_t = obj.simulate_1run( ad , M1 , 'in_x0' , in_x0  )
			%	
			%
			%Inputs:
			%	ad - An Aff_Dyn object.
			%		 This should be the affine dynamics for which this controller was synthesized.
			%		 If it is not, then no guarantees on performance can be maintained.
			%	M1 - A nonnegative scalar.
			%	in_str - This input string should tell what type of constant you would like to hold.
			%
			%Outputs:
			%	x_0_t - An n x (T+1) matrix which defines the trajectory of the state
			%			for a feasible realization of the random variables associated
			%			with the tuple (ad,M1,sigma)

			%++++++++++++++++
			%Input Processing
			obj = varargin{1};
			ad = varargin{2};
			M1 = varargin{3};

			if nargin < 4
				rand_word_ind = randi(length(obj.L),1);
				sig = obj.L{rand_word_ind};
				T = length(obj.L{rand_word_ind});
			else
				if isa(varargin{4},'char')
					%This means that the newer version of simulate_1run is being called.
					%Read the characters.
					flag_ind = 4;
					while flag_ind <= nargin
						switch varargin{flag_ind}
						case 'in_sigma'
							if (size(varargin{flag_ind+1},1) ~= 1) %|| (size(in_sig,2) ~= T)
								error('Input word is not a single word or does not have the correct length.' )
							end
							sig = varargin{flag_ind+1};
							T = length(sig);
							%Increment Flag Index
							flag_ind = flag_ind+2;
						case 'in_w'
                            %Choose a random word to make our T
                            rand_word_ind = randi(length(obj.L),1);
                            sig = obj.L{rand_word_ind};
                            T = length(sig);
							if ~exist('T')
								warning('simulate_1run was called with ''in_w'' parameter flag, but ''T'' was not defined.')
							end
							in_w = varargin{flag_ind+1};
							if ((size(in_w,1) ~= size(ad.B_w,2)) || (size(in_w,2) ~= T))
								error('The dimensions of the input w sequence are not correct.')
							end
							%Increment Flag Index
							flag_ind = flag_ind+2;
						case 'in_x0'
							x0 = varargin{flag_ind+1};
							if ((size(x0,1) ~= size(ad.A,2)) || (size(x0,2) ~= 1))
								error('The dimensions of the input w sequence are not correct.')
							end
							%Increment Flag Index
							flag_ind = flag_ind+2;
						otherwise
							error(['Unrecognized input to simulate_1run: ' varargin{flag_ind} ])
						end
					end
				elseif isa(varargin{4},'double')
					%This means that the older version of simulate_1run is being called.
					%Read the double.
					in_sig = varargin{4};
				else
					error('Unrecognized type for fourth argument.')
				end
			end

			%+++++++++
			%Algorithm

			%Constants
			% T = size(obj.L,2);
			n = size(ad.A,1);
			wd = size(ad.B_w,2);
			vd = size(ad.C_v,2);

			%Generate Random Variables
			if ~exist('x0')
				x0 = unifrnd(-M1,M1,n,1);
			end

			if exist('in_w')
				w = in_w;
			else
				if ~isnan(ad.eta_w)
					w  = unifrnd(-ad.eta_w,ad.eta_w,wd,T);
				elseif isa(ad.P_w,'Polyhedron')
					w = ad.P_w.V'*unifrnd(0,1,size(ad.P_w.V,1),T);
				else
					error('Do not understand how to handle this definition of w.')
				end
			end

			if ~isnan(ad.eta_v)
				v  = unifrnd(-ad.eta_v,ad.eta_v,vd,T);
			elseif isa(ad.P_v,'Polyhedron')
				v = ad.P_v.V'*unifrnd(0,1,size(ad.P_v.V,1),T)
			else
				error('Do not understand how to handle this definition of v.')
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

		function [run_data_x,run_data_x_norm] = simulate_n_runs( varargin )
			%Description:
			%	Uses the information provided in the affine dynamics instance ad and the problem
			%	specification M1 along with the given controller to simulate as many runs as the user
			%	would like.
			%
			%Usage:
			%	x_tensor = simulate_n_runs( obj , ad , M1 , num_runs , in_sig )
			%
			%Inputs:
			%	
			%
			%Outputs:
			%	run_data_x -	A 1 x 'num_runs' cell array of n x (T_i+1) matrices which defines
			%					the num_runs trajectories of the state.
			%			   		T_i can change for each run.
			%	run_data_x_norm - 	A cell array of

			%++++++++++++++++
			%Input Processing

			if nargin < 4
				error('Not enough inputs.')
			end

			obj = varargin{1};
			ad = varargin{2};
			M1 = varargin{3};
			num_runs = varargin{4};

			if nargin >= 5
				if isa(varargin{5},'char')
					%This means that the newer version of simulate_1run is being called.
					%Read the characters.
					flag_ind = 5;
					while flag_ind <= nargin
						switch varargin{flag_ind}
						case 'in_sigma'
							if (size(varargin{flag_ind+1},1) ~= 1) %|| (size(in_sig,2) ~= T)
								error('Input word is not a single word or does not have the correct length.' )
							end
							in_sig = varargin{flag_ind+1};
							T = length(in_sig);
							%Increment Flag Index
							flag_ind = flag_ind+2;
						case 'in_w'
							in_w = varargin{flag_ind+1};
							if ((size(in_w,1) ~= size(ad.B_w,2)) || (size(in_w,2) ~= T))
								error('The dimensions of the input w sequence are not correct.')
							end
							%Increment Flag Index
							flag_ind = flag_ind+2;
						otherwise
							error(['Unrecognized input to simulate_n_runs: ' varargin{flag_ind} ])
						end
					end
				elseif isa(varargin{5},'double')
					%This means that the older version of simulate_1run is being called.
					%Read the double.
					in_sig = varargin{5};
				else
					error('Unrecognized type for fourth argument.')
				end
			end

			%+++++++++
			%Algorithm

			%% Constants
			n = size(ad.A,1);

			%% Implement Loop
			run_data_x = {};
			run_data_x_norm = {};
			for run_ind = 1:num_runs
				%
				if nargin < 5
					run_data_x{run_ind} = obj.simulate_1run( ad , M1 );
				elseif exist('in_sig')
					run_data_x{run_ind} = obj.simulate_1run( ad , M1 , 'in_sigma', in_sig );
				elseif exist('in_w')
					run_data_x{run_ind} = obj.simulate_1run( ad , M1 , 'in_sigma' , in_sig , 'in_w' , in_w )
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
