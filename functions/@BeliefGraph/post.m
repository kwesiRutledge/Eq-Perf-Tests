function ancest_nodes = post(varargin)
	%Description:
	%	Identifies what nodes could possibly arise after reaching the current Belief Node according to the dynamics
	%	given in lcsas.
	%
	%Usage:
	%	lcsas.post(BN,P_u,P_x0)
	%	lcsas.post(BN,P_u,P_x0,'debug',debug_flag)
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

	lcsas 	= varargin{1};
	BN 		= varargin{2};
	P_u 	= varargin{3};
	P_x0 	= varargin{4};

	if ~isa(lcsas,'LCSAS')
		error('You are using a deprecated function call. Please make sure that the input to post is an LCSAS object.')
	end

	if ~isa(BN,'BeliefNode')
		error('You are expected to provide a belief node as the second input.')
	end

	if nargin > 4
		arg_idx = 5;
		while arg_idx <= nargin
			switch varargin{arg_idx}
				case 'debug'
					debug_flag = varargin{arg_idx+1};
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

	%Get All Combinations of the node's subset
	node_p_set = BN.idx_powerset_of_subL();

	n = size(lcsas.Dyn(1).A,1);
	m = size(lcsas.Dyn(1).B,2);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Define Consistency Sets %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	Phi_sets = {}; visible_transitions = [1:length(node_p_set)];
	for p_set_idx = 1:length(node_p_set)
		Phi_sets{p_set_idx} = lcsas.consistent_set(BN.t+1,{BN.subL.words{node_p_set{p_set_idx}}},P_u,P_x0);
		%If any of the Phi's are empty,
		%then it is impossible for a transition to exist between the node c_level(node_ind) and the node associated with Phi
		if Phi_sets{p_set_idx}.isEmptySet
			visible_transitions(p_set_idx) = 0;
		end
	end 

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Identify Which Consistency Sets Can Be Independently Detected %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%For each possible transition, see if its transition set is completely contained by another transition set
	for ut_idx = 1:length(node_p_set)
		if visible_transitions(ut_idx) ~= 0
			Projx_Phi1 = [eye(n*((BN.t+1)+1)+m*(BN.t+1)) zeros(n*((BN.t+1)+1)+m*(BN.t+1),Phi_sets{ut_idx}.Dim-n*((BN.t+1)+1)-m*(BN.t+1))]*Phi_sets{ut_idx};
			%vt_tilde = visible_transitions(visible_transitions ~= ind_ut);
			for ch_idx = [ut_idx+1:length(node_p_set)]
				if debug_flag >= 1
					disp(['ch_idx = ' num2str(ch_idx)])
				end
				Projx_Phi2 = [eye(n*((BN.t+1)+1)+m*(BN.t+1)) zeros(n*((BN.t+1)+1)+m*(BN.t+1),Phi_sets{ch_idx}.Dim-n*((BN.t+1)+1)-m*(BN.t+1))]*Phi_sets{ch_idx};
				temp_diff = Projx_Phi1 \ Projx_Phi2;
				%Check to see if temp_diff is a single Polyhedron/Polytope (if it is multiple then the set difference is not empty)
				if length(temp_diff) == 1
					if (temp_diff.isEmptySet) && (~(Projx_Phi1.isEmptySet)) && (~(Projx_Phi2.isEmptySet))
						if debug_flag >= 1
							disp(['temp_diff.isEmptySet = ' num2str(temp_diff.isEmptySet) ' for:'])
							disp(['- Phi1(' num2str([BN.subL.words{node_p_set{ut_idx}}]) ')' ])
							disp(['- Phi2(' num2str([BN.subL.words{node_p_set{ch_idx}}]) ')' ])
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
		temp_L = Language();
		temp_L.words = {BN.subL.words{node_p_set{visible_transitions(trans_idx)}}};

		c_node = BeliefNode(temp_L,BN.t+1);
		ancest_nodes = [ancest_nodes,c_node];
	end

end