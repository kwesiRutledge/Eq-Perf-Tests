classdef BeliefNode
	%Description:
	%
	%Properties:
	%	- subL:	A Language object that represents the belief at the current time step.
	%			i.e. One of the mode sequences in this language is occurring, if we have
	%			reached this node.
	%
	%	- t:	Time at which the belief is held.
	%			Assumed to be between 0 and the length of the longest word in subset_L.
	%
	%	- c_set:	A Polyhedron that represents what trajectories (y_[0:t],u_[0:t-1])
	%	- FullTrajectorySet: 
	%
	%Methods:
	%	- idx_powerset_of_subL
	%	- post
	%	- is_eq
	%	- find_longest_horizon
	%

	properties(SetAccess = protected)
		c_set_present;
		phi_set_present;
	end

	properties
		subL;
		t;
		c_set;
		FullTrajectorySet; %Also known as Phi
	end

	methods
		function BN = BeliefNode(varargin)
			%Description:
			%	A Belief Node is defined to contain:
			%		- A Belief on which mode sequences are being executed (represented as a subset of the
			%		  language L defined for the LCSAS), and
			%		- A time at which that belief is held
			%Usage:
			%	BN = BeliefNode(subset_L,t0)
			%	BN = BeliefNode(subset_L,t0,c_set)
			%	BN = BeliefNode(subset_L,t0,c_set,phi_set)
			%	BN = BeliefNode(subset_L,t0,'FullTrajectorySet',phi_set)

			debug_flag = 0;

			%%%%%%%%%%%%%%%%%%%%%%
			%% Input Processing %%
			%%%%%%%%%%%%%%%%%%%%%%

			if nargin < 2
				error('Not enough input arguments.')
			end

			subset_L = varargin{1};
			t0 = varargin{2};

			varargin_idx = 3;
			while varargin_idx < nargin
				if isa(varargin{varargin_idx},'Polyhedron')
					switch varargin_idx
						case 3
							c_set = varargin{varargin_idx};
						case 4
							phi_set = varargin{varargin_idx};
						otherwise
							error('Unrecognized BeliefNode Constructor.')
					end
					continue;
				end


				switch varargin{varargin_idx}
					case 'FullTrajectorySet'
						phi_set = varargin{varargin_idx+1};
						varargin_idx = varargin_idx+2;
					otherwise
						error(['Unrecognized BeliefNode Constructor input: ' varargin{varargin_idx} ])
				end

			end

			% switch nargin
			% 	case 2
			% 		%warning('This method for initializing Belief Nodes is deprecated. Please use the 3 argument version.')
			% 		1;
			% 	case 3
			% 		c_set = varargin{3};
			% 	case 4
			% 		c_set = varargin{3};
			% 		phi_set = varargin{4};
			% 	otherwise
			% 		error(['It is not allowable to call this class with ' num2str(nargin) ' arguments.' ])
			% end

			if ~(iscell(subset_L) || isa(subset_L,'Language') )
				error('Expected subset of language to be a cell array or a Language object.')
			end

			if iscell(subset_L)
				warning('The version of this class that uses L as a cell array is deprecated. We cannot guarantee that functions will work anymore as the class has been modified to support Language objects.')
				%Convert subset_L to a Language
				temp_arr = subset_L;
				clear subset_L
				subset_L = Language();
				subset_L.words = temp_arr;
			end

			%%%%%%%%%%%%%%%
			%% Constants %%
			%%%%%%%%%%%%%%%

			if ~exist('c_set')
				c_set = [];
			end

			if ~exist('phi_set')
				phi_set = [];
			end

			%%%%%%%%%%%%%%%%%%%%%%%%%
			%% Assigning Variables %%
			%%%%%%%%%%%%%%%%%%%%%%%%%

			BN.subL = subset_L;
			BN.t = t0;

			BN.c_set = c_set;
			BN.FullTrajectorySet = phi_set;

			BN.c_set_present = ~isempty(c_set);
			BN.phi_set_present = ~isempty(phi_set);

		end

		function subsets = idx_powerset_of_subL(obj)
			%Description:
			%	Assigns to each word in subL.words an index.
			%	Then returns all possible subsets of the INDICES.
			%
			
			node_p_set = {};
			for comb_length = 1:length(obj.subL.words)
				temp_combs = nchoosek([1:length(obj.subL.words)],comb_length);
				for comb_ind = 1:size(temp_combs,1)
					node_p_set{end+1} = temp_combs(comb_ind,:);
				end
			end
			subsets = node_p_set;
		end

		function ancest_nodes = post(varargin)
			%Description:
			%	Identifies what nodes could possibly arise after reaching the current Belief Node according to the dynamics
			%	given in lcsas.
			%
			%Usage:
			%	post(BN,lcsas,P_u,P_x0)
			%	post(BN,lcsas,P_u,P_x0,'debug',debug_flag)
			%	post(BN,lcsas,P_u,P_x0,'fb_method','state')
			%
			%Inputs:
			%	lcsas - An array of Aff_Dyn() objects.
			%			 

			%%%%%%%%%%%%%%%%%%%%%%
			%% Input Processing %%
			%%%%%%%%%%%%%%%%%%%%%%

			if nargin < 4
				error('Not enough input arguments.')
			end

			BN 		= varargin{1};
			lcsas 	= varargin{2};
			P_u 	= varargin{3};
			P_x0 	= varargin{4};

			if ~isa(lcsas,'LCSAS')
				error('You are using a deprecated function call. Please make sure that the input to post is an LCSAS object.')
			end

			if nargin > 4
				arg_idx = 5;
				while arg_idx <= nargin
					switch varargin{arg_idx}
						case 'debug'
							debug_flag = varargin{arg_idx+1};
							arg_idx = arg_idx + 2;
						case 'fb_method'
							fb_method = varargin{arg_idx+1};
							arg_idx = arg_idx + 2;
						otherwise
							error(['Unexpected input: ' varargin{arg_idx}])
					end
				end
			end

			%%%%%%%%%%%%%%%
			%% Constants %%
			%%%%%%%%%%%%%%%

			if ~exist('debug_flag')
				debug_flag = 0;
			end

			if ~exist('fb_method')
				fb_method = 'output';
			end

			%Get All Combinations of the node's subset
			%node_p_set = BN.idx_powerset_of_subL();
			subL_p_set = BN.subL.powerset();

			n = size(lcsas.Dyn(1).A,1);
			m = size(lcsas.Dyn(1).B,2);

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%% Define Consistency Sets %%
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			Phi_sets = {};
			Consist_sets = {};
			visible_transitions = [1:length(subL_p_set)];
			for p_set_idx = 1:length(subL_p_set)
				disp(['p_set_idx = ' num2str(p_set_idx)])
				%Check first to see if any of the words in this sublanguage are too short to create something at time t+1
				if subL_p_set(p_set_idx).find_shortest_length() >= BN.t + 1
					[ Consist_sets{p_set_idx} , ~ ] = lcsas.consistent_set(BN.t+1,subL_p_set(p_set_idx),P_u,P_x0,'fb_method',fb_method);
					%If any of the Phi's are empty,
					%then it is impossible for a transition to exist between the node c_level(node_ind) and the node associated with Phi
					if Consist_sets{p_set_idx}.isEmptySet
						visible_transitions(p_set_idx) = 0;
					end
				else
					visible_transitions(p_set_idx) = 0;
				end
					
			end 

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%% Identify Which Consistency Sets Can Be Independently Detected %%
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			%For each possible transition, see if its transition set is completely contained by another transition set
			for ut_idx = 1:length(subL_p_set)
				if debug_flag >= 1
					disp(['ut_idx = ' num2str(ut_idx)])
				end
				if visible_transitions(ut_idx) ~= 0
					for ch_idx = [ut_idx+1:length(subL_p_set)]
						if debug_flag >= 1
							disp(['ch_idx = ' num2str(ch_idx)])
						end
						temp_diff = Consist_sets{ut_idx} \ Consist_sets{ch_idx};
						%Check to see if temp_diff is a single Polyhedron/Polytope (if it is multiple then the set difference is not empty)
						if length(temp_diff) == 1
							if (temp_diff.isEmptySet) && (~(Consist_sets{ut_idx}.isEmptySet)) && (~(Consist_sets{ch_idx}.isEmptySet))
								if debug_flag >= 2
									disp(['temp_diff.isEmptySet = ' num2str(temp_diff.isEmptySet) ' for:'])
									disp(['- Phi1(' num2str([subL_p_set(ut_idx).words{:}]) ')' ])
									disp(['- Phi2(' num2str([subL_p_set(ch_idx).words{:}]) ')' ])
									disp(' ')
								end
								visible_transitions(ut_idx) = 0;
							% elseif (ch_idx == ind_ut) && (~Projx_Phi1.isEmptySet)
							% 	visible_transitions(ind_ut) = ch_idx;
							end 
						end
					end
				end
			end
			%Process visible_transitions matrix
			visible_transitions = visible_transitions(visible_transitions ~= 0);
			visible_transitions = unique(visible_transitions);

			%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%% Create Ancestor Nodes %%
			%%%%%%%%%%%%%%%%%%%%%%%%%%%
			ancest_nodes = [];
			for trans_idx = 1:length(visible_transitions)
				temp_L = subL_p_set(visible_transitions(trans_idx));
				temp_consist_set = Consist_sets{visible_transitions(trans_idx)};

				c_node = BeliefNode(temp_L,BN.t+1,temp_consist_set);
				ancest_nodes = [ancest_nodes,c_node];
			end

		end

		function eq_flag = is_eq(obj,BN)
			%Description:
			%	Returns true if this belief node (obj) is equal to the belief node BN.

			if obj.subL.is_eq(BN.subL) && (obj.t == BN.t)
				eq_flag = true;
			else
				eq_flag = false;
			end

		end

		function longest_T = find_longest_horizon(obj)
			%Description:
			%	Searches through all elements of the subL.words for this node and determines
			%	how much longer of a future the system can have (at maximum).

			longest_T = -1;

			for L_idx = 1:length(obj.subL.words)
				longest_T = max(length(obj.subL.words{L_idx}),longest_T) - obj.t;
			end
		end
	end

end