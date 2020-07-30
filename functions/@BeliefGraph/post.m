function ancest_nodes = post(varargin)
	%Description:
	%	Identifies what nodes could possibly arise after reaching the current Belief Node according to the dynamics
	%	given in lcsas.
	%
	%Usage:
	%	ancest_nodes = BG.post(BN,P_u,P_x0)
	%	ancest_nodes = BG.post(BN,P_u,P_x0,'debug',debug_flag)
	%	
	%Assumption:
	%	This function assumes that the BeliefGraph function contains the following member variables
	%		- UsedProjection: This flag indicates whether or not the method will use projection or not.
	%		- FeedbackMethod: Indicates if the system is using state or output feedback.
	%		- UsedAcceleratedAlgorithms: Indicates if the accelerated "post" algorithms should be used.
	%		- UsedUnobservabilityChecks: Indicates that the "post" algorithm checks for unobservable edges/transitions.
	%
	%Inputs:
	%	BG 			- A Belief Graph object.
	%	lcsas 		- An array of Aff_Dyn() objects.
	%	P_u 		- A Polyhedron() object that defines the input constraints.
	%				  The input at each time must lie within the polyhedron.
	%	P_x0 		- A Polyhedron() object that defines the set of states from which the initial state is contained.
	%	debug_flag  - A nonnegative integer indicating the level of debug information that the user wants.

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if nargin < 4
		error('Not enough input arguments.')
	end

	BG 		= varargin{1};
	BN 		= varargin{2};
	P_u 	= varargin{3};
	P_x0 	= varargin{4};

	if ~isa(BN,'BeliefNode')
		error('Please make sure that the second input to post is a BeliefNode object.')
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

	%%%%%%%%%%%%%%
	%% Defaults %%
	%%%%%%%%%%%%%%

	if ~exist('debug_flag')
		debug_flag = 0;
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	lcsas = BG.lcsas;
	subL = BN.subL;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Select which Post function to call %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if BG.UsedProjection
		ancest_nodes = BG.post_proj(BN,P_u,P_x0, ...
									'debug',debug_flag, ...
									'fb_method',BG.FeedbackMethod , ...
									'accel_flag' , BG.UsedAcceleratedAlgorithms , ...
									'use_unobservability_checks' , BG.UsedUnobservabilityChecks );
		return
	else
		ancest_nodes = BG.post_noproj(BN,P_u,P_x0, ...
										'debug',debug_flag, ...
										'fb_method',BG.FeedbackMethod , ...
										'accel_flag' , BG.UsedAcceleratedAlgorithms);
		return
	end

end