classdef FHAE_pb
%Description:
%	This class is the prefix-based controller that was developed in my most recent write up. (#31)
	properties
		L;
		F_set;
		u0_set;
		BG;
	end

	methods
		%Constructor
		function contr = FHAE_pb( varargin )
			%Description:
			%	Simply copy everything over
			%
			%Usage:
			%	contr = FHAE_pb( L  , F_set , u0_set )
			%	contr = FHAE_pb( BG , F_set , u0_set )
			%
			%Inputs:
			%	L - 	A cell array of one dimension.
			%	F_set - A 1 x |L| cell array of feedback matrices which are indexed
			%			based on the index of the matching word in L.
			
			%%%%%%%%%%%%%%%%%%%%%%
			%% Input Processing %%
			%%%%%%%%%%%%%%%%%%%%%%

			if isa(varargin{1},'BeliefGraph')
				BG = varargin{1};
				F_set = varargin{2};
				u0_set = varargin{3};
			else
				L = varargin{1};
				F_set = varargin{2};
				u0_set = varargin{3};

				if isnumeric(L)
					%Convert L to cell array.
					L_temp = {};
					for word_ind = 1:size(L,1)
						L_temp{word_ind} = L(word_ind,:);
					end
					contr.L = Language(L_temp);
				elseif iscell(L) 
					contr.L = Language(L);
				elseif isa(L,'Language')
					contr.L = L;
				else
					error('Unknown type of L.')
				end
			end

			%Get any additional inputs.
			varargin_idx = 4;
			while(varargin_idx <= nargin)
				switch varargin{varargin_idx}
					case 'fb_type'
						error('This flag is currently not in use.')
						fb_type = varargin{varargin_idx};
						if ~(strcmp(fb_type,'disturbance') || strcmp(fb_type,'output'))
							error('Unexpected feedback type. Expecting ''disturbance'' or ''output''.')
						end
						varargin_idx = varargin_idx + 2;
					otherwise
						error(['Unrecognized input string: ' varargin{varargin_idx}])
				end
			end

			%%%%%%%%%%%%%%%
			%% Constants %%
			%%%%%%%%%%%%%%%

			%%%%%%%%%%%%%%%
			%% Algorithm %%
			%%%%%%%%%%%%%%%

			if exist('L')
				%If the first input was a lnaguage, check the value of L.
				contr.L = L;
				contr.BG = [];
			elseif exist('BG')
				contr.L = [];
				contr.BG = BG;
			end
				
			contr.F_set = F_set;
			contr.u0_set = u0_set;
		end

		%Apply the correct control
		function u = apply_control( varargin )
			%FHAE_pb.apply_control
			%Description:
			%	Receives the currently measured output values
			%	and computes the controller's output value.
			%
			%Usage:
			%	u = apply_control(obj , observed_w , y_vec)
			%	u = apply_control(obj , y_vec)

			%%%%%%%%%%%%%%%%%%%%%%
			%% Input Processing %%
			%%%%%%%%%%%%%%%%%%%%%%

			if (nargin < 2) || (nargin > 3)
				error('Expected 2 or 3 inputs.')
			end

			obj = varargin{1};
			switch nargin
				case 2
					y_vec = varargin{2};
					mode_prop = 'unobserved';
				case 3
					observed_w = varargin{2};
					y_vec = varargin{3};
					mode_prop = 'observed';
				otherwise
					error('Expected 2 or 3 inputs.')
			end

			%%%%%%%%%%%%%%%%%%%%%%
			%% Useful Constants %%
			%%%%%%%%%%%%%%%%%%%%%%

			ow_len = length(observed_w);
			num_words = length(obj.L);
			m = size(obj.F_set{1},1) / length(obj.L{1});
			p = size(obj.F_set{1},2) / length(obj.L{1});

			%%%%%%%%%%%%%%%
			%% Algorithm %%
			%%%%%%%%%%%%%%%

			switch mode_prop
				case 'observed'
					%Match observed_w with a word from L
					temp_L = Language(obj.L);
					first_match_ind = temp_L.find_a_word_with_pref(observed_w);

					%Obtain the correct feedback matrices
					F_t  = obj.F_set{first_match_ind}([m*(ow_len-1)+1:m*ow_len],[1:p*ow_len]);
					u0_t = obj.u0_set{first_match_ind}([m*(ow_len-1)+1:m*ow_len]);
				otherwise
					body
			end

			% Compute Output
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
