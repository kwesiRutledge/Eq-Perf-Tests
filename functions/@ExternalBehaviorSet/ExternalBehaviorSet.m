classdef ExternalBehaviorSet < handle
	%Description:
	%	The external behavior set associated with a given LCSAS and knowledge sequence.
	%	Also known as a Consistency Set.
	%
	%To-do:
	%	- Add second constructor which does what get_closed_loop_consistent_internal_behavior_set_matrices.m is doing now.

	properties
		System;
		AsPolyhedron;
		t;
		KnowledgeSequence;

		% Parent Set
		ParentInternalBehaviorSet;

		% Other Things
		Dim;
		Settings;
	end

	methods

		function [ebs] = ExternalBehaviorSet(varargin)
			%Description:
			%	Defines an object which represents a polyhedron of external behaviors that are consistent for modes of a
			%	LCSAS.
			%	To place the system into clearer focus. We have a Language-Constrained Switched Affine System (LCSAS):
			%
			%	x_{t+1} = A_{q_t} x_t + B_{q_t} u_t + f_{q_t} + w_t
			%	
			%	where:
			%			- q_t is a natural number that describes the current mode at time t
			%			- w_t belongs to the set W_{q_t} which varies with the mode
			%
			%	Under output feedback, the consistency set can also be written as:
			%					{ [y]  | }
			%		C(\sigma) =	{ [u]  | }
			%	
			%	Under state feedback, the consistency set can also be written as:
			%					{ [x]  | }
			%		C(\sigma) =	{ [u]  | }
			%
			%Inputs:
			%	lcsas 		- An array of Aff_Dyn() objects. Hopefully the dimensions are all appropriately checked so that
			%				  the state is the proper size in all Aff_Dyn(). Inputs as well, etc.
			%	t 			- The time of interest
			%	L 			- The set of words under consideration.
			%				  We would like to find the set of states for which it is possible to reach when under ALL switching
			%				  sequences in this set with the same inputs (albeit with different disturbances).
			%	use_proj 	- Boolean (true or false).
			%				  Used to tell the function to either skip the creation of Consist_set (false) or complete the
			%				  computation of Consist_set (true) which requires projection operations to be called and may be very slow.
			%
			%Example Usage:
			%	[ebs] = ExternalBehaviorSet(lcsas,KnowledgeSequence)
			%	ebs = ExternalBehaviorSet(lcsas,KnowledgeSequence,'fb_method','state')			

			%% Input Processing
			[ lcsas , KnowledgeSequence , ebs_settings ] = input_processing_ExternalBehaviorSet(varargin{:});
            ebs.Settings = ebs_settings;

            %% Constants %%
            [ n_x , n_u , n_y , n_w , n_v ] = lcsas.Dimensions();
            
			%% Create Internal Behavior Set
			ebs.ParentInternalBehaviorSet = InternalBehaviorSet(lcsas, KnowledgeSequence,'ebs_settings_struct',ebs_settings);
            t = ebs.ParentInternalBehaviorSet.t;
            
			%% Get The Matrice

			%% Define Object
			ebs.t = t;
			ebs.System = lcsas;
            ebs.KnowledgeSequence = KnowledgeSequence;
			ebs.AsPolyhedron = [];
            
            switch ebs_settings.fb_type
                case 'state'
                    ebs.Dim = n_x * (t+1) + n_u * t;
                case 'output'
                    ebs.Dim = n_y * (t+1) + n_u * t;
                otherwise
                    error(['Unexpected feedback type given in ExternalBehaviorSet(): ' ebs_settings.fb_type ])
            end
                

		end

		function [poly_out] = ToPolyhedron(ebs)
			%Description:
			%	Computes the MPT3 Polyhedron() representation of the ExternalBehaviorSet.

			%% Constants
			System = ebs.System;
			t = ebs.t;
			ebs.ParentInternalBehaviorSet;

			L = System.L;

			[ n_x , n_u , n_y , n_w , n_v ] = System.Dimensions();

			%% Algorithm
			if ~isempty(ebs.AsPolyhedron)
				%Polyhedron version has already been computed.
				poly_out = ebs.AsPolyhedron;
			else
				
				%Check to see if ibs has computed Polyhedron representation.
				if isempty(ebs.ParentInternalBehaviorSet.AsPolyhedron)
					ebs.ParentInternalBehaviorSet.ToPolyhedron();
				end

				poly_out = [ eye(n_x*(t+1) + n_u*t), zeros(n_x*(t+1) + n_u*t, ebs.ParentInternalBehaviorSet.Dim - n_x*(t+1) - n_u*t ) ] * ebs.ParentInternalBehaviorSet.AsPolyhedron;

				ebs.AsPolyhedron = poly_out;

			end



		end

		function [tf] = le(ebs1,ebs2)
			%Description:
			%	Defines the <= operator.
			%	This will define set inclusion of the left (ebs1) by the right (ebs2) sets. It transforms each set
			%	into a Polyhedron and then does the comparison.

			%% Constants

			%% Algorithm
			if isempty(ebs1.AsPolyhedron)
				poly1 = ebs1.ToPolyhedron();
			else
				poly1 = ebs1.AsPolyhedron;
			end 

			if isempty(ebs2.AsPolyhedron)
				poly2 = ebs2.ToPolyhedron();
			else
				poly2 = ebs2.AsPolyhedron;
			end

			tf = (poly1 <= poly2); %Use MPT3's inclusion operator.

		end

		function [tf] = ge(ebs1,ebs2)
			%Description:
			%	Defines the >= operator.
			%	This will define the set inclusion of the right (ebs2) by the left (ebs1). It transforms each set into
			%	a Polyhedron and then does the comparison.

			tf = le(ebs2,ebs1);

		end

		function [tf_matrix] = find_coverage_relationship(ebs_array)
			%Description:
			%	Determines which of the ExternalBehaviorSet objects in ebs_array
			%	covers which other items using a simple containment check.
			%
			%Usage:
			%	tf_matrix = ebs_array.find_coverage_relationship()

			%% Constants

			num_ebs = length(ebs_array);

			%% Algorithm
			
			%Fill In Diagonal
			tf_matrix = false(num_ebs);
			for ebs_index1 = 1:num_ebs
				tf_matrix(ebs_index1,ebs_index1) = true;
			end

			all_ebs_indices = [1:num_ebs];

			for ebs_index1 = 2:num_ebs
				for ebs_index2 = find(all_ebs_indices ~= ebs_index1)
					% Identify if ebs(index1) is contained by ebs(index2)
					tf_matrix( ebs_index1 , ebs_index2 ) = (ebs_array(ebs_index1) <= ebs_array(ebs_index2)); 
				end
			end

		end

		function [ empty_flags ] = IsEmpty( ebs )
			%Description:
			%	This function identifies which of the ExternalBehaviorSets are empty.
			%	If a given path's internal behavior set is empty, then we will mark that set appropriately.

			%% Input Checking

			% if isscalar(ebs)
			% 	error('The matrices A,b, Ae, and be are expected to be numeric in order to check to see if this is an empty set.')
			% end

			%% Variables

            empty_flags = logical.empty;
            
			%% Algorithm

			if isscalar(ebs)

				empty_flags = ebs.ParentInternalBehaviorSet.IsEmpty();

			elseif isvector(ebs) && iscell(ebs)
                
				for ebs_index = 1:length(ebs)
					temp_single_ebs = ebs{ebs_index};
					empty_flags = [ empty_flags ; temp_single_ebs.IsEmpty() ];
				end

			elseif isvector(ebs) && ~iscell(ebs)

				for ebs_index = 1:length(ebs)
					temp_single_ebs = ebs(ebs_index);
					empty_flags = [ empty_flags ; temp_single_ebs.IsEmpty() ];
				end

			else
				error('Input to IsEmpty() must be a scalar or vector.')
			end

		end

		function [ variables_out , constraints_out ] = CreateContainmentConstraint( ebs_in , ebs_circum )
			%Description:
			%	Creates a containment constraint between the two external behavior sets ebs_in (the in-body, or the 
			%	set that we hope to contain) and ebs_circum (the circumbody, or the set that we hope contains the other).

			% Constants
			cg = constr_gen(0); %Create constraint generator.

			Dim = ebs_in.Dim;

			ibs_in = ebs_in.ParentInternalBehaviorSet;
			ibs_circum = ebs_circum.ParentInternalBehaviorSet;

			% Algorithm

			A_in = [ ibs_in.A ; ibs_in.Ae ; -ibs_in.Ae ];
			b_in = [ ibs_in.b ; ibs_in.be ; -ibs_in.be ];

			A_circum = [ ibs_circum.A ; ibs_circum.Ae ; -ibs_circum.Ae ];
			b_circum = [ ibs_circum.b ; ibs_circum.be ; -ibs_circum.be ];

			[variables_out, constraints_out] = cg.create_sadraddini_AH_inclusion_constr( ...
													zeros(Dim,1) , ibs_in.SelectExternalBehavior() , 	A_in , b_in , ...
													zeros(Dim,1) , ibs_circum.SelectExternalBehavior() , A_circum , b_circum );

		end

	end

end

function [ lcsas , KnowledgeSequence , cs_settings ] = input_processing_ExternalBehaviorSet(varargin)
	%Description:
	%	Processes the inputs given to ExternalBehaviorSet() constructor.

	% Check for Minimum Number of Arguments.
	if nargin < 2
		error('Not enough input arguments.')
	end

	% Parse first 5 Arguments
	lcsas = varargin{1};
	KnowledgeSequence = varargin{2};
	
	% Check if Some Arguments are of Appropriate Type

	if ~isa(lcsas,'LCSAS')
		error('Expecting the first input to be a LCSAS object.')
	end
	lcsas.check('U','X0','L')

	%% Default Values

	cs_settings = struct( ...
		'fb_type', 'state', ...
		'use_proj', true, ...
		'reduce_flag', true, ...
		'OpenLoopOrClosedLoop', 'Open', ...
		'K', [], ...
		'k', [] ...
		);

	varargin_idx = 3;
	while varargin_idx <= nargin
		switch varargin{varargin_idx}
			case {'fb_method','fb_type'}
				cs_settings.fb_type = varargin{varargin_idx+1};
				if ~(strcmp(cs_settings.fb_type,'state') || strcmp(cs_settings.fb_type,'output'))
					error(['Invalid feedback type: ' fb_type ])
				end
				varargin_idx = varargin_idx + 2;
			case 'use_proj'
				cs_settings.use_proj = varargin{varargin_idx+1};
				if ~islogical( cs_settings.use_proj )
					error('The flag for ''use_proj'' should be a boolean.')
				end
				varargin_idx = varargin_idx + 2;
			case 'reduce_representation'
				cs_settings.reduce_flag = varargin{varargin_idx+1};
				if ~islogical( cs_settings.reduce_flag )
					error('The flag for ''reduce_flag'' should be a boolean.')
				end
				varargin_idx = varargin_idx + 2;
			case 'OpenLoopOrClosedLoop'
				cs_settings.OpenLoopOrClosedLoop = varargin{varargin_idx+1};
				if strcmp(cs_settings.OpenLoopOrClosedLoop,'Closed')
					cs_settings.K = varargin{varargin_idx+2};
					cs_settings.k = varargin{varargin_idx+3};
					varargin_idx = varargin_idx + 4;
				else
					varargin_idx = varargin_idx + 2;
				end
			otherwise
				error(['Unexpected additional input:' varargin{varargin_idx}])
		end
	end

	if (length(KnowledgeSequence) < 0)
		error(['KnowledgeSequence must have a length greater than 0.'])
	end


end