classdef ConsistentBeliefsController < handle
%Description:
%	
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
		K_set;	%Linear Gains
		k_set;	%Offset Gain

		System;

		% Knowledge Sequences Replace the Need for a Belief Graph
		KnowledgeSequences;
		ConsistencySets; 	%Each set should correspond to the consistency set
		
		% External Behavior history
		u_hist;
		x_hist;
		
		y_hist;

		b_hist;

		% Settings
		Settings;
	end

	methods
		%Constructor
		function cbc = ConsistentBeliefsController( varargin )
			%Description:
			%	This class is the prefix-based controller that uses Consistency Sets according to the
			%	method in our LCSS work.
			%
			%Usage:
			%	contr = ConsistentBeliefsController( lcsas , BeliefSequences , K_set , k_set )
			%
			%Inputs:
			%	BG - 	A BeliefGraph() object representing how the belief might change over time.
			%	F_set - A 1 x |L| cell array of feedback matrices which are indexed
			%			based on the index of the matching word in L.
			
			%%%%%%%%%%%%%%%%%%
			%% Set Defaults %%
			%%%%%%%%%%%%%%%%%%

			cbc.Settings = struct('FeedbackMethod','Disturbance (State)');

			%%%%%%%%%%%%%%%%%%%%%%
			%% Input Processing %%
			%%%%%%%%%%%%%%%%%%%%%%

			if (~isa(varargin{1},'LCSAS'))
				error('The input system must be an LCSAS object.')
			end

			cbc.System = 				varargin{1};
			cbc.KnowledgeSequences = 	varargin{2};
			cbc.K_set = 				varargin{3};
			cbc.k_set = 				varargin{4};

			%Get any additional inputs.
			argument_index = 5;
			while(argument_index <= nargin)
				switch varargin{argument_index}
					case 'FeedbackType'
						cbc.Settings.FeedbackMethod = varargin{argument_index};
						argument_index = argument_index + 2;
					otherwise
						error(['Unrecognized input string: ' varargin{argument_index} ' to ConsistentBeliefsController().'])
				end
			end

			%%%%%%%%%%%%%%%
			%% Constants %%
			%%%%%%%%%%%%%%%

			num_sequences = size(cbc.KnowledgeSequences,2);
			TimeHorizon = length(cbc.System.L.words{1});

			%%%%%%%%%%%%%%%
			%% Algorithm %%
			%%%%%%%%%%%%%%%

			ConsistencySets = {};

			for t = 0:TimeHorizon-1
				for ks_index = 1:num_sequences
					ConsistencySets{t+1,ks_index} = ExternalBehaviorSet( cbc.System , cbc.KnowledgeSequences([1:t+1],ks_index), ...
                                                                        'OpenLoopOrClosedLoop','Closed',cbc.K_set{ks_index}, cbc.k_set{ks_index} );
				end
			end
			cbc.ConsistencySets = ConsistencySets;

			cbc.u_hist = []; cbc.x_hist = []; cbc.y_hist = [];
			cbc.b_hist = [];
		end

		function eb = history_to_external_behavior( cbc )
			%Description:
			%	Stack the values which determine external behavior depending on the settings of the controller.

			switch cbc.Settings.FeedbackMethod
			case 'Disturbance (State)'
				eb = [ 	cbc.history_to_x_vec() ; 
						reshape(cbc.u_hist, prod(size(cbc.u_hist)) , 1 ) ];
			otherwise
				error(['Unexpected feedback method given to create_external_behavior.'])
			end
		end

		function x_vec = history_to_x_vec(cbc)
			%Description:
			%	Stacks the values of the x_hist vector into a single vector according in chronological order.
			%
			%Usage:
			%	x_vec = cbc.history_to_x_vec();

			x_vec = reshape( cbc.x_hist , prod( size(cbc.x_hist) ) , 1 );

		end

		function y_vec = history_to_y_vec(cbc)
			%Description:
			%	Stacks the values of the y_hist vector into a single vector according in chronological order.
			%
			%Usage:
			%	y_vec = cbc.history_to_y_vec();

			y_vec = reshape( cbc.y_hist , prod( size(cbc.y_hist) ) , 1 );

		end

		function w_vec = history_to_w_vec(cbc)
			%Description:
			%	Identifies the sequence of disturbances w_vec 
			%
			%Usage:
			%	w_vec = cbc.history_to_w_vec()

			%% Constants %%

			System = cbc.System;

			t = size(cbc.u_hist,2);
			% most_recent_hypotheses = cbc.b_hist(end);

			%% Algorithm %%

			if t == 0
				error(['history_to_w_vec() is not meant to be called with t=0.'])
			end

			w_vec = [];

			for tau = 0:t-1
				x_tau 	= cbc.x_hist(:,tau+1);
				x_taup1 = cbc.x_hist(:,tau+1+1);
				u_tau 	= cbc.u_hist(:,tau+1);
                
                hypothesis_at_t = cbc.b_hist(tau+1+1);

				first_hypothesis = hypothesis_at_t.words{1};
				hypothesisDynamicsAtTau = System.Dyn( first_hypothesis(tau+1) );

				w_tau = hypothesisDynamicsAtTau.find_w_that_completes( x_tau , u_tau , x_taup1 );
				w_vec = [w_vec; w_tau];

			end
		end

		function clear_histories(cbc)
			%Description:
			%	Clear all the history variables.
			%
			%Usage:
			%	cbc.clear_histories()

			cbc.u_hist = [];
			cbc.x_hist = [];
			
			cbc.y_hist = [];

			cbc.b_hist = [];
		end

		%Apply the correct control
		function u = compute_control( cbc )
			%FHAE_pb.compute_control
			%Description:
			%	Receives the currently measured output values
			%	and computes the controller's output value.
			%
			%Usage:
			%	u = compute_control( cbc )
			%	u = cbc.compute_control()
			
			%%%%%%%%%%%%%%%%%%%%%%
			%% Useful Constants %%
			%%%%%%%%%%%%%%%%%%%%%%

			System = cbc.System;
			Settings = cbc.Settings;
			[ n_x , n_u , n_y , n_w , n_v ] = System.Dimensions();
			% p = size(y_mat,1);
			% m = size(obj.F_set{1},1) / T1;

			KnowledgeSequences = cbc.KnowledgeSequences;
			[T,num_sequences] = size(KnowledgeSequences);
			t = size(cbc.u_hist,2);

			%%%%%%%%%%%%%%%
			%% Algorithm %%
			%%%%%%%%%%%%%%%

			switch Settings.FeedbackMethod
			case 'Disturbance (State)'

				% Prefix Detection
				gain_idx = cbc.prefix_detection();
                cbc.b_hist = [cbc.b_hist; KnowledgeSequences(t+1,gain_idx)];

				if t == 0

					k_t = cbc.k_set{gain_idx}(n_u*t+[1:n_u]);
					
					u = k_t;
				else	

					w_vec = cbc.history_to_w_vec();

					%Obtain the correct feedback matrices
					K_t = cbc.K_set{gain_idx}([n_u*t+1:n_u*(t+1)],[1:n_w*t]);
					k_t = cbc.k_set{gain_idx}([n_u*t+1:n_u*(t+1)]);

					% Compute Output
					u = K_t * w_vec + k_t;
				end

			otherwise
				error(['Unexpected feedback method (' Settings.FeedbackMethod ') given to compute_control().' ])
			end

			cbc.u_hist = [cbc.u_hist,u];
		end

		function [ detected_index ] = prefix_detection( cbc )
			%Description:
			%	Given the histories observed so far, find the path/KnowledgeSequence which matches our
			%	currently observed external behavior the best.

			%% Constants %%

			KnowledgeSequences = cbc.KnowledgeSequences;
			u_hist = cbc.u_hist;
			b_hist = cbc.b_hist;

			eb0 = cbc.history_to_external_behavior();
			t = size(u_hist,2);
            
			num_sequences = size(KnowledgeSequences,2);

			%% Algorithm %%
			candidate_indices = [];

			%Only consider indices which have the right prefix.
			if ~isempty(b_hist)
				for knowl_seq_index = 1:num_sequences
					temp_prefix = KnowledgeSequences([1:t+1-1],knowl_seq_index);
					if temp_prefix == cbc.b_hist
						candidate_indices = [candidate_indices, knowl_seq_index];
					end
				end
			else
				candidate_indices = [1:num_sequences];
			end

			matching_indices = [];
			for knowl_seq_index = candidate_indices
				temp_cs = cbc.ConsistencySets{t+1,knowl_seq_index};
				if temp_cs.contains(eb0)
					matching_indices = [matching_indices; knowl_seq_index];
				end
            end

            if length(matching_indices) == 0
                disp('There was an issue identifying the external behavior!')
                disp(['eb0 = ' num2str(eb0')])
                if t > 0
                    disp(['cbc.b_hist ='])
                    for tau = 0:length(cbc.b_hist)-1
                        disp(cbc.b_hist(tau+1))
                    end
                end
            end
            
			%Search through all matching indices for the one with maximum cardinality.
			detected_index = matching_indices(1);
			detected_prefix = KnowledgeSequences([1:t+1],detected_index);

			for mi_index = 2:length(matching_indices)
				mi = matching_indices(mi_index);
				temp_prefix = KnowledgeSequences([1:t+1],mi);
				if temp_prefix >= detected_prefix
					detected_prefix = temp_prefix;
					detected_index = mi;
				end
			end


		end

		%Simulation of a single run
		function [ x_0_t, u_0_tm1 , y_0_t , sig ] = simulate_1run( cbc )
			%Description:
			%	Uses the information provided in the Affine Dynamics instance ad and
			%	the problem specification M1 along with the given controller
			%	to simulate a single run of the prefix based controller.
			%Usages:
			%	x_0_t = obj.simulate_1run()
			%
			%Inputs:
			%
			%Outputs:
			%	x_0_t - An n x (T+1) matrix which defines the trajectory of the state
			%			for a feasible realization of the random variables associated
			%			with the tuple (ad,M1,sigma)

			%%%%%%%%%%%%%%%
			%% Constants %%
			%%%%%%%%%%%%%%%
			
			System = cbc.System;
			L = System.L;
			[ n_x , n_u , n_y , n_w , n_v ] = System.Dimensions();

			x_0_t = []; y_0_t = []; sig = []; u_0_tm1 = [];
			in_w = []; x0 = [];

			if nargin < 2
				rand_word_ind = randi(L.cardinality(),1);
				sig = L.words{rand_word_ind};
				T = length(sig);
			else
				[sig,T,in_w,x0] = cbc.simulate_1run_input_helper( varargin );
			end

			%%%%%%%%%%%%%%%
			%% Algorithm %%
			%%%%%%%%%%%%%%%
            
			%Constants
			% T = size(obj.L,2);
            
			%Generate Random Variables
			x0 = sample_once_from( System.X0 );

			if ~isempty(in_w)
				w = in_w;
			else
				w = [];
				for symbol_index = 1:length(sig)
					temp_symbol = sig(symbol_index);
					W_si = System.Dyn( temp_symbol ).P_w;

					w = [ w , sample_once_from(W_si) ];
				end
            end
            
			v = [];
			for symbol_index = 1:length(sig)
				temp_symbol = sig(symbol_index);
				V_si = System.Dyn( temp_symbol ).P_v;

				v = [ v , sample_once_from(V_si) ];
			end

			% v = obj.gen_rand_vars( lcsas , sig , 'v' );

			%Initilize system
			x_t = x0;
			q_t = sig(1);

			y_t = System.Dyn(q_t).C*x0 + System.Dyn(q_t).C_v*v(:,1);

            %Initialize Controller Histories
            cbc.clear_histories();
			x_0_t = x_t; cbc.x_hist = x_0_t;
			y_0_t = y_t; cbc.y_hist = y_0_t;
			x_tp1 = Inf(n_x,1);
            
			for t = 0:T-1
				q_t = sig(t+1);

				%Use Affine Dynamics with proper control law.
				x_tp1 = System.Dyn(q_t).A * x_t + ...
						System.Dyn(q_t).B * cbc.compute_control() + ...
						System.Dyn(q_t).B_w * w(:,t+1) + ...
						System.Dyn(q_t).f ;

				%Update other variables in system
				x_t = x_tp1;
				x_0_t = [x_0_t, x_t];
				cbc.x_hist = x_0_t;

				if( t == T-1 )
					continue;
				end
				q_tp1 = sig(t+2);

				y_t = System.Dyn(q_tp1).C*x_t + System.Dyn(q_tp1).C_v*v(:,t+1);
				y_0_t = [y_0_t, y_t];
				cbc.y_hist = y_0_t;
            end

            if T == 3
                disp('What?')
            end
            
			u_0_tm1 = cbc.u_hist;

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
			%	x_tensor = simulate_n_runs( obj , M1 , num_runs , in_sig )
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

			error('This function has not been rewritten for LCSAS.')

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

		function plotConsistencySets( cbc )
			%Description:
			%	Plots the consistency sets in a single figure.

			%% Constants %%

			[ T , num_sequences ] = size(cbc.ConsistencySets);
			System = cbc.System;

			[ n_x , n_u , n_y , n_w , n_v ] = System.Dimensions();

			%% Algorithm %%

			% Plotting each of the slices for each time t
			for sequence_index = 1:num_sequences
				figure;
				title(['Consistency Sets for Sequence #' num2str(sequence_index) '/' num2str(num_sequences) ])
				for t = 0:T-1
					subplot(1,T,t+1)
					plot(cbc.ConsistencySets{t+1,sequence_index}.ToPolyhedron().projection(n_x*t+[1:n_x]))
					title(['t=' num2str(t) ])
				end
				
			end

		end

		function [ plotHandles , R ] = plotReachableSets( varargin )
			%Description:
			%	Tries to plot each of the reachable sets of the controller.
			%Usage:
			%	cbc.plotReachableSets()
			%	cbc.plotReachableSets('PlotTarget',true,'TargetSet',X_T)

			%% Constants %%

			cbc = varargin{1};

			[ T , num_sequences ] = size(cbc.ConsistencySets);
			System = cbc.System;
			ConsistencySets = cbc.ConsistencySets;

			[ n_x , n_u , n_y , n_w , n_v ] = System.Dimensions();
			L = System.L;

			x0 = System.X0.V';

			%% Input Processing %%

			pRS_settings = struct(	'PlotOnCurrentFigure',false, ...
									'PlotTarget',false, ...
									'TargetSet',[], ...
									'TargetColor', 'yellow' , ...
									'TargetAlpha', 0.5 );

			tempColorArray = {};
			for sequence_index = 1:num_sequences
				tempColorArray{sequence_index} = 'cyan';
			end
			pRS_settings.ReachableSetColors = tempColorArray;

			if nargin > 1
				argument_index = 2;
				while nargin >= argument_index
					switch varargin{argument_index}
					case 'ReachableSetColors'
						pRS_settings.ReachableSetColors = varargin{argument_index+1};

						if length(pRS_settings.ReachableSetColors) ~= num_sequences
							error(['Expected for there to be enough colors in ReachableSetColors for each sequence in cbc.KnowledgeSequences. There are ' num2str(length(ReachableSetColors)) 'colors and ' num2str(num_sequences) 'sequences.'  ])
						end

						argument_index = argument_index + 2;
					case 'PlotTarget'
						pRS_settings.PlotTarget = varargin{argument_index+1};
						argument_index = argument_index + 2;
					case 'TargetSet'
						pRS_settings.TargetSet = varargin{argument_index+1};
						argument_index = argument_index + 2;
					otherwise
						error(['Unexpected input to plotReachableSets: ' varargin{argument_index} ])
					end
				end
			end

			%% Algorithm %%

			% Create Word Dependent Disturbance sets

			WT_i = {};
			for sequence_index = [1:num_sequences]
				% Collect Consistency Sets
				temp_sequence = cbc.KnowledgeSequences(:,sequence_index);
				ibs_ts = InternalBehaviorSet(System,temp_sequence,'OpenLoopOrClosedLoop','Closed',cbc.K_set{sequence_index},cbc.k_set{sequence_index});
				WT_i_part = ibs_ts.SelectW()*ibs_ts.ToPolyhedron();

				% Append an Extra W term to the end of WT_i_part.
				L_Tm1 = temp_sequence(end);
				word1 = L_Tm1.words{1};
				[~,word1_index] = L.contains(word1);

				WT_i{sequence_index} = WT_i_part * System.Dyn(word1_index).P_w;

			end

			% Create Reachable Sets Using Affine Expressions
			R = {};

			for sequence_index = [1:num_sequences]
				temp_sequence = cbc.KnowledgeSequences(:,sequence_index);
				L_Tm1 = temp_sequence(end);

				[S_w_i,S_u_i,~,J_i,f_bar_i] = get_mpc_matrices(System,'word',L_Tm1.words{1});

				K = cbc.K_set{sequence_index}; k = cbc.k_set{sequence_index};
				K_trimmed = K([1:n_u*T],[1:n_w*T]);
				k_trimmed = k([1:n_u*T],1);


				S_w_trimmed = S_w_i([1:n_x*(T+1)],[1:n_w*T]);
				S_u_trimmed = S_u_i([1:n_x*(T+1)],[1:n_u*T]);
				J_trimmed = J_i([1:n_x*(T+1)],:);
				f_bar_prime = S_w_trimmed*f_bar_i([1:n_w*T]);

				R{sequence_index} = ( S_w_trimmed + S_u_trimmed * K_trimmed)*WT_i{sequence_index} + S_u_trimmed*k_trimmed + J_trimmed * x0 + f_bar_prime ;
			end

			% ========
			% Plotting
			% ========

			plotHandles = {};

			for sequence_index = 1:num_sequences

				pH_i = [];

				figure;
				hold on;

				% Plot with target if asked
				if pRS_settings.PlotTarget

					if isempty(pRS_settings.TargetSet)
						error(['Target set was not provided to the function plotReachableSets.'])
					end

					pH_i = plot(pRS_settings.TargetSet, 'color', pRS_settings.TargetColor , 'alpha' , pRS_settings.TargetAlpha);

				end

				pH_i(end+1) = scatter(x0(1),x0(2));

				for t = 1:T
					% Plot The Polyhedron For Word 1
					if length(plotHandles) < 3
						pH_i(end+1) = plot( R{sequence_index}.projection(n_x*t+[1:n_x]), ...
							'color',pRS_settings.ReachableSetColors{sequence_index});
					else
						plot( R{sequence_index}.projection(n_x*t+[1:n_x]), ...
							'color',pRS_settings.ReachableSetColors{sequence_index} );
					end
				end

				title(['Reachable Sets for Mode Sequence #' num2str(sequence_index) ])

				%Create legend
				if pRS_settings.PlotTarget
					legend(pH_i,'$$\mathcal{X}_T$$','$$x_0$$',['$$R(\mathbf{m}^{(' num2str(sequence_index)' ')})$$'],'Interpreter','latex');
				else
					legend(pH_i,'$$x_0$$',['$$R(\mathbf{m}^{(' num2str(sequence_index)' ')})'],'Interpreter','latex');
				end

				plotHandles{sequence_index} = pH_i;

			end

		end

	end

end
