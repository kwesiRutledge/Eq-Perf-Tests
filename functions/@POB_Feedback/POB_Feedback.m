classdef POB_Feedback < handle
%Description:
%	This class is the prefix-based, belief space feedback that was developed in the paper
%	'Belief-prefix Control for Autonomously Dodging Switching Disturbances' (#60)
%
%Member Functions:
%	apply_control
%	simulate_1run
%	gen_rand_vars
%	simulate_1run_input_helper
%	simulate_n_runs
%
%Usage:
%
	properties
		F_set;
		u0_set;

		BG;
		u_hist;
		b_hist;
	end

	methods
		%Constructor
		function contr = POB_Feedback( varargin )
			%Description:
			%	This class is the prefix-based, belief space feedback that was developed in the paper
			%	'Belief-prefix Control for Autonomously Dodging Switching Disturbances' (#60)
			%
			%Usage:
			%	contr = POB_Feedback( BG , F_set , u0_set )
			%
			%Inputs:
			%	BG - 	A BeliefGraph() object representing how the belief might change over time.
			%	F_set - A 1 x |L| cell array of feedback matrices which are indexed
			%			based on the index of the matching word in L.
			
			%%%%%%%%%%%%%%%%%%%%%%
			%% Input Processing %%
			%%%%%%%%%%%%%%%%%%%%%%

			BG = varargin{1};
			F_set = varargin{2};
			u0_set = varargin{3};

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

			contr.BG = BG;

			contr.F_set = F_set;
			contr.u0_set = u0_set;

			contr.u_hist = [];
			contr.b_hist = [];
		end

		%Apply the correct control
		function u = apply_control( varargin )
			%FHAE_pb.apply_control
			%Description:
			%	Receives the currently measured output values
			%	and computes the controller's output value.
			%
			%Usage:
			%	u = apply_control(obj , observed_w , y_mat)
			%	u = apply_control(obj , y_mat)
			

			%%%%%%%%%%%%%%%%%%%%%%
			%% Input Processing %%
			%%%%%%%%%%%%%%%%%%%%%%

			if (nargin < 2) || (nargin > 2)
				error('Expected 2 inputs.')
			end

			obj = varargin{1};
			y_mat = varargin{2};
			mode_prop = 'unobserved';

			%%%%%%%%%%%%%%%%%%%%%%
			%% Useful Constants %%
			%%%%%%%%%%%%%%%%%%%%%%

			p = size(y_mat,1);

			num_words = length(obj.BG.BeliefLanguage.words);
			T1 = length(obj.BG.BeliefLanguage.words{1});
			m = size(obj.F_set{1},1) / T1;
			%p = size(obj.F_set{1},2) / T1;

			y_vec = reshape(y_mat,prod(size(y_mat)),1);
			num_obsvs = length(y_vec)/p;

			%%%%%%%%%%%%%%%
			%% Algorithm %%
			%%%%%%%%%%%%%%%

			gain_idx = obj.prefix_detection( y_mat );
			obj.BG.BeliefLanguage.words{gain_idx};

			%Obtain the correct feedback matrices
			F_t  = obj.F_set{gain_idx}([m*(num_obsvs-1)+1:m*num_obsvs],[1:p*num_obsvs]);
			u0_t = obj.u0_set{gain_idx}([m*(num_obsvs-1)+1:m*num_obsvs]);

			% Compute Output
			u = F_t * y_vec + u0_t;

			obj.u_hist = [obj.u_hist,u];
		end

		%Simulation of a single run
		function [ x_0_t, u_0_tm1 , y_0_t , sig ] = simulate_1run( varargin )
			%Description:
			%	Uses the information provided in the Affine Dynamics instance ad and
			%	the problem specification M1 along with the given controller
			%	to simulate a single run of the prefix based controller.
			%Usages:
			%	x_0_t = obj.simulate_1run( ad , M1 , sig )
			%	x_0_t = obj.simulate_1run( ad , M1 , 'in_sigma' , sig )
			%	x_0_t = obj.simulate_1run( ad , M1 , 'in_sigma' , sig , 'in_w' , in_w )
			%	x_0_t = obj.simulate_1run( ad , M1 , 'in_x0' , in_x0  )
			%	x_0_t
			%
			%Inputs:
			%	ad - An Aff_Dyn object.
			%		 This should be the affine dynamics for which this controller was synthesized.
			%		 If it is not, then no guarantees on performance can be maintained.
			%	lcsas - An LCSAS object.
			%
			%	M1 - A nonnegative scalar.
			%	in_str - This input string should tell what type of constant you would like to hold.
			%
			%Outputs:
			%	x_0_t - An n x (T+1) matrix which defines the trajectory of the state
			%			for a feasible realization of the random variables associated
			%			with the tuple (ad,M1,sigma)

			%%%%%%%%%%%%%%%%%%%%%%
			%% Input Processing %%
			%%%%%%%%%%%%%%%%%%%%%%
			
			obj = varargin{1};
			in_sys = varargin{2};
			M1 = varargin{3};

			if (~isa(in_sys,'LCSAS'))
				error('The input system must be an LCSAS object.')
			end

			x_0_t = []; y_0_t=[]; sig = []; u_0_tm1 = [];

			%Save system in a new container
			lcsas = in_sys;

			in_w = []; x0 = [];

			if nargin < 4
				rand_word_ind = randi(length(lcsas.L.words),1);
				sig = lcsas.L.words{rand_word_ind};
				T = length(sig);
			else
				[sig,T,in_w,x0] = obj.simulate_1run_input_helper( varargin );
			end

			%%%%%%%%%%%%%%%
			%% Algorithm %%
			%%%%%%%%%%%%%%%

			%Constants
			% T = size(obj.L,2);
			n = size(lcsas.Dyn(1).A,1);
			wd = size(lcsas.Dyn(1).B_w,2);
			vd = size(lcsas.Dyn(1).C_v,2);

			%Generate Random Variables
			if isempty(x0)
				x0 = obj.gen_rand_vars( lcsas , sig , 'x0' , M1 );
			end

			if ~isempty(in_w)
				w = in_w;
			else
				w = obj.gen_rand_vars( lcsas , sig , 'w' );
			end

			v = obj.gen_rand_vars( lcsas , sig , 'v' );

			%Simulate system forward.
			x_t = x0;
			q_t = sig(1);

			y_t = lcsas.Dyn(q_t).C*x0 + lcsas.Dyn(q_t).C_v*v(:,1);
			
			x_0_t = x_t;
			y_0_t = y_t;
			x_tp1 = Inf(n,1);

			for t = 0:T-1
				q_t = sig(t+1);

				%Use Affine Dynamics with proper control law.
				x_tp1 = lcsas.Dyn(q_t).A * x_t + ...
						lcsas.Dyn(q_t).B * obj.apply_control( y_0_t ) + ...
						lcsas.Dyn(q_t).B_w * w(:,t+1) + ...
						lcsas.Dyn(q_t).f ;
				%Update other variables in system
				x_t = x_tp1;
				x_0_t = [x_0_t, x_t];

				if( t == T-1 )
					continue;
				end
				q_tp1 = sig(t+2);

				y_t = lcsas.Dyn(q_tp1).C*x_t + lcsas.Dyn(q_tp1).C_v*v(:,t+1);
				y_0_t = [y_0_t, y_t];
			end

			u_0_tm1 = obj.u_hist;

		end

		%Generate random variables
		function [rand_var] = gen_rand_vars( varargin )
			%Description:
			%	Uses the dynamics and the mode signal to generate a sequence of disturbance variables.
			%
			%Usage:
			%	w = FHAE_pb.gen_rand_vars( aff_dyn , word , 'w' )
			%	v = FHAE_pb.gen_rand_vars( aff_dyn , word , 'v' )
			%	x0 = FHAE_pb.gen_rand_vars( aff_dyn , word , 'x0' , M1 )
			%	x0 = FHAE_pb.gen_rand_vars( aff_dyn , word , 'x0' , Px0 )
			%	
			%	w = FHAE_pb.gen_rand_vars( lcsas , word ,'w' )
			%	w = FHAE_pb.gen_rand_vars( lcsas , word ,'v' )
			%	w = FHAE_pb.gen_rand_vars( lcsas , word ,'x0' , M1 )


			%%%%%%%%%%%%%%%%%%%%%%
			%% Input Processing %%
			%%%%%%%%%%%%%%%%%%%%%%

			obj = varargin{1};
			in_sys = varargin{2};
			word = varargin{3};
			disturb_flag = varargin{4};

			%%%%%%%%%%%%%%%
			%% Constants %%
			%%%%%%%%%%%%%%%

			rand_var = NaN;
			T = length(word);

			%%%%%%%%%%%%%%%%
			%% Algorithms %%
			%%%%%%%%%%%%%%%%

			obj = varargin{1};
			in_sys = varargin{2};

			if isa(in_sys,'Aff_Dyn')

				%% Constants %%
				ad = in_sys;

				n = size(ad.A,1);
				wd = size(ad.B_w,2);
				vd = size(ad.C_v,2);
				
				%% Algorithm %%
				switch disturb_flag
				case 'w'
					if ~isnan(in_sys.eta_w)
						rand_var  = unifrnd(-in_sys.eta_w,in_sys.eta_w,wd,T);
					elseif isa(in_sys.P_w,'Polyhedron')
						num_verts = size(lcsas.Dyn(1).P_w.V,2);
						convex_comb = unifrnd(0,1,num_verts,T);
						convex_comb = convex_comb./repmat(sum(convex_comb),num_verts,1);

						rand_var = in_sys.P_w.V'*convex_comb;
					else
						error('Do not understand how to handle this definition of w.')
					end
				case 'v'
					%Generate the measurement disturbances.
					if ~isnan(ad.eta_v)
						rand_var  = unifrnd(-ad.eta_v,ad.eta_v,vd,T);
					elseif isa(ad.P_v,'Polyhedron')
						%Create Convex Combination Vector
						num_verts = size(lcsas.Dyn(1).P_v.V,1);
						convex_comb = unifrnd(0,1,num_verts,T);
						convex_comb = convex_comb./repmat(sum(convex_comb),num_verts,1);

						rand_var = ad.P_v.V'*convex_comb;
					else
						error('Do not understand how to handle this definition of v.')
					end
				case 'x0'
					x0_set_def = varargin{5};
					%Generate x0
					if isa(x0_set_def,'Polyhedron')
						%Input Processing
						Px0 = x0_set_def;

						%Convex Combination
						convex_comb = unifrnd(0,1,size(Px0.V,1),T);
						convex_comb = convex_comb./sum(convex_comb);

						rand_var = Px0.V'*convex_comb; %'
					elseif isscalar(x0_set_def)
						M1 = x0_set_def;
						rand_var = unifrnd(-M1,M1,n,1);
					else
						error('Do not understand how to handle this definition of x0.')
					end
				otherwise
					error('Unexpected disturbance flag.')
				end

			elseif isa(in_sys,'LCSAS')
				
				%% Constants %%
				lcsas = in_sys;

				n = size(lcsas.Dyn(1).A,1);
				wd = size(lcsas.Dyn(1).B_w,2);
				vd = size(lcsas.Dyn(1).C_v,2);
				
				%% Algorithm %%
				switch disturb_flag
				case 'w'
					if ~isnan(lcsas.Dyn(1).eta_w)
						rand_var = [];
						for symb_idx = 1:T
							symb_t = word(symb_idx);
							w_t = unifrnd(	-lcsas.Dyn(symb_t).eta_w, ...
											lcsas.Dyn(symb_t).eta_w, wd,1);

							rand_var  = [rand_var, w_t];
						end
					elseif isa(lcsas.Dyn(1).P_w,'Polyhedron')
						rand_var = [];

						convex_comb = unifrnd(0,1,size(lcsas.Dyn(1).P_w.V,1),T);
						num_verts = size(lcsas.Dyn(1).P_w.V,1);
						convex_comb = convex_comb./repmat(sum(convex_comb),num_verts,1);

						for symb_idx = 1:T
							symb_t = word(symb_idx);
							w_t = lcsas.Dyn(symb_t).P_w.V'*convex_comb(:,symb_idx); %'

							rand_var = [rand_var, w_t];
						end
					else
						error('Do not understand how to handle this definition of w.')
					end
				case 'v'
					%Generate the measurement disturbances.
					if ~isnan(lcsas.Dyn(1).eta_v)
						rand_var = [];
						for symb_idx = 1:T
							symb_t = word(symb_idx);
							v_t = unifrnd(	-lcsas.Dyn(symb_t).eta_v, ...
											lcsas.Dyn(symb_t).eta_v, vd,1);

							rand_var  = [rand_var, v_t];
						end
					elseif isa(lcsas.Dyn(1).P_v,'Polyhedron')
						rand_var = [];

						convex_comb = unifrnd(0,1,size(lcsas.Dyn(1).P_v.V,1),T);
						num_verts = size(lcsas.Dyn(1).P_v.V,1);
						convex_comb = convex_comb./repmat(sum(convex_comb),num_verts,1);

						for symb_idx = 1:T
							symb_t = word(symb_idx);
							v_t = lcsas.Dyn(symb_t).P_v.V'*convex_comb(:,symb_idx); %'

							rand_var = [rand_var, v_t];
						end
					else
						error('Do not understand how to handle this definition of v.')
					end
				case 'x0'
					x0_set_def = varargin{5};
					%Generate x0
					if isa(x0_set_def,'Polyhedron')
						%Input Processing
						Px0 = x0_set_def;

						%Convex Combination
						convex_comb = unifrnd(0,1,size(Px0.V,1),T);
						convex_comb = convex_comb/sum(convex_comb);

						rand_var = Px0.V'*convex_comb; %'
					elseif isscalar(x0_set_def)
						M1 = x0_set_def;
						rand_var = unifrnd(-M1,M1,n,1);
					else
						error('Do not understand how to handle this definition of x0.')
					end
				otherwise
					error('Unexpected disturbance flag.')
				end

			else
				error('Input system can only be an Aff_Dyn or LCSAS object.')
			end

		end

		function [sig,T,in_w,x0] = simulate_1run_input_helper( obj, arg_arr )
			%Description:
			%
			%Usage:
			%	[sig,T,~,~] = FHAE_pb.simulate_1run_input_helper( arg_arr )
			%	[sig,T,in_w,~] = FHAE_pb.simulate_1run_input_helper( arg_arr )
			%	[sig,T,in_w,x0] = FHAE_pb.simulate_1run_input_helper( arg_arr )

			%%%%%%%%%%%%%%%%%%%%%%
			%% Input Processing %%
			%%%%%%%%%%%%%%%%%%%%%%

			obj = arg_arr{1};
			in_sys = arg_arr{2};
			M1 = arg_arr{3};

			%%%%%%%%%%%%%%%
			%% Algorithm %%
			%%%%%%%%%%%%%%%
			sig = []; T = []; in_w = []; x0 = [];

			flag_ind = 4;

			if isa(in_sys,'Aff_Dyn')

				ad = in_sys;

				while flag_ind <= length(arg_arr)
					switch arg_arr{flag_ind}
					case 'in_sigma'
						if (size(arg_arr{flag_ind+1},1) ~= 1) %|| (size(in_sig,2) ~= T)
							error('Input word is not a single word or does not have the correct length.' )
						end
						sig = arg_arr{flag_ind+1};
						if ~exist('T')
							T = length(in_sig);
						end
						%Increment Flag Index
						flag_ind = flag_ind+2;
					case 'in_w'
	                    in_w = arg_arr{flag_ind+1};
						if (size(in_w,1) ~= size(ad.B_w,2))
							error('The dimensions of the input w sequence are not correct.')
						end
						if ~exist('T')
							T = size(in_w,2);
						end
						%Increment Flag Index
						flag_ind = flag_ind+2;
					case 'in_x0'
						x0 = arg_arr{flag_ind+1};
						if ((size(x0,1) ~= size(ad.A,2)) || (size(x0,2) ~= 1))
							error('The dimensions of the input w sequence are not correct.')
						end
						%Increment Flag Index
						flag_ind = flag_ind+2;
					otherwise
						error(['Unrecognized input to simulate_1run: ' arg_arr{flag_ind} ])
					end
				end
			elseif isa(in_sys,'LCSAS')

				lcsas = in_sys;

				while flag_ind <= length(arg_arr)
					switch arg_arr{flag_ind}
					case 'in_sigma'
						if (size(arg_arr{flag_ind+1},1) ~= 1) %|| (size(in_sig,2) ~= T)
							error('Input word is not a single word or does not have the correct length.' )
						end
						sig = arg_arr{flag_ind+1};
						T = length(sig);
						%Increment Flag Index
						flag_ind = flag_ind+2;
					case 'in_w'
						in_w = arg_arr{flag_ind+1};
						if (size(in_w,1) ~= size(lcsas.Dyn(1).B_w,2)) %|| (size(in_w,2) ~= T))
							error('The dimensions of the input w sequence are not correct.')
						end

						if ~exist('T')
							T = size(in_w,2);
						end
						%Increment Flag Index
						flag_ind = flag_ind+2;
					case 'in_x0'
						x0 = arg_arr{flag_ind+1};
						if ((size(x0,1) ~= size(lcsas.Dyn(1).A,2)) || (size(x0,2) ~= 1))
							error('The dimensions of the input x0 value is not correct.')
						end
						%Increment Flag Index
						flag_ind = flag_ind+2;
					otherwise
						error(['Unrecognized input to simulate_1run: ' arg_arr{flag_ind} ])
					end
				end
			else
				error('in_sys is expected to be either an Aff_Dyn or LCSAS object.')
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
							if ~exist('T')
								T = length(in_sig);
							end
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
