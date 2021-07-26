classdef InternalBehaviorSet < handle
	%Description:
	%	The internal behavior set associated with a given LCSAS.
	%
	%To-do:
	%	- Add second constructor which does what get_closed_loop_consistent_internal_behavior_set_matrices.m is doing now.
	%	- Add empty initializer for IBS and only create one when ToPolyhedron() is called.

	properties
		System;
		AsPolyhedron;
		t;
		KnowledgeSequence;

		% H Representation
		A;
		b;
		Ae;
		be;

		% Associated Gains
		K;
		k;

		% Other Things
		Dim;
		ibs_settings;
	end

	methods
		function [ibs] = InternalBehaviorSet(varargin)
			%internal_behavior_set.m
			%Description:
			%	Finds a polyhedron that describes what pairs of state, input, and disturbance sequences are compatible/feasible from ALL
			%	switching sequences defined by the sequence KnowledgeSequence.
			%
			%	The KnowledgeSequence is a vector describing the following sequence of estimates about the true state of the system;
			%
			%		KnowledgeSequence = [ L_{0} , L_{1} , ... , L_{t-1} , L_t ]
			%	
			%	where L_{i} indicates the estimate of the mode at the time i.
			%	This means that length(KnowledgeSequence) == t.
			%
			%	To place the system into clearer focus. We have a Language-Constrained Switched Affine System (LCSAS):
			%
			%	x_{t+1} = A_{q_t} x_t + B_{q_t} u_t + f_{q_t} + w_t
			%	
			%	where:
			%			- q_t is a natural number that describes the current mode at time t
			%			- w_t belongs to the set W_{q_t} which varies with the mode
			%	The consistency set can also be written as
			%
			%									  { [x]  | }
			%		IB_lcsas(KnowledgeSequence) = { [u]  | }
			%									  { [w]  | }
			%									  { [x0] | }
			%		
			%										{ [y]  | }
			%										{ [u]  | }
			%		IB_lcsas(KnowledgeSequence) = 	{ [w]  | }
			%										{ [v]  | }
			%										{ [x0] | }
			%										{ [x]  | }
			%
			%Inputs:
			%	lcsas 		- An array of Aff_Dyn() objects. Hopefully the dimensions are all appropriately checked so that
			%				  the state is the proper size in all Aff_Dyn(). Inputs as well, etc.
			%	use_proj 	- Boolean (true or false).
			%				  Used to tell the function to either skip the creation of Consist_set (false) or complete the
			%				  computation of Consist_set (true) which requires projection operations to be called and may be very slow.
			%
			%Example Usage:
			%	[ internal_behavior_set ] = InternalBehaviorSet( lcsas , KnowledgeSequence )
			%	[ internal_behavior_set ] = InternalBehaviorSet( lcsas , KnowledgeSequence , 'fb_type' , 'state' )
			%
			%Assumptions:
			%	This formulation assumes that the system does not include a disturbed measurements. i.e. We can perfectly observe the state
			%	at each time.

			%% Input Processing
			[ lcsas , KnowledgeSequence , K , k , ibs_settings ] = ibs.ip_InternalBehaviorSet( varargin{:} );
			
			ibs.System = lcsas;
			ibs.KnowledgeSequence = KnowledgeSequence;
			ibs.K = K;
			ibs.k = k;
			ibs.ibs_settings = ibs_settings;

			%% Optional Method to Exit The Constructor Early
			if ibs_settings.ReturnEarly
				ibs.A = ibs_settings.A; ibs.b = ibs_settings.b;
				ibs.Ae = ibs_settings.Ae; ibs.be = ibs_settings.be;

				ibs.Dim = ibs_settings.Dim;

				ibs.System = lcsas;
				ibs.KnowledgeSequence = KnowledgeSequence;

				return;
			end

			%% Constants

			KS_len = length(KnowledgeSequence);
			t = KS_len-1;
			LargestLangInSequence = KnowledgeSequence(1);
			T = length(lcsas.L.words{1});

		    [ n_x , n_u , n_y , n_w , n_v ] = lcsas.Dimensions();
		    X0 = lcsas.X0;
		    P_u = lcsas.U;

			%% Algorithm

			if T <= t
				error(['The length of the word is ' num2str(T) ' and yet the knowledge sequence would like to know about time ' num2str(t) ' which is impossible. Make sure that T > t.' ])
			end

			ibs.t = t;
		        
			% Debugging Strategy
			if ibs_settings.debug > 0
				disp('Debugging info...')
			end

			A = []; b = [];
			Ae = []; be = [];

		    for tau = 1:t
				L_tau = KnowledgeSequence(tau+1);

				%Create a Part of the A Matrix for Each word in L_tau
				[A_projected,b_projected,Ae_projected,be_projected] = ibs.CreatePolytopeMatricesAtTime(tau,L_tau,ibs_settings);
				%Extend all projected matrices so that they can fit in A,b,Ae,be
				[ tempA , tempb , tempAe , tempbe ] = ibs.AdjustIBSMatricesFor( tau , L_tau , ...
																				A_projected , b_projected , Ae_projected , be_projected );

				A = [ A ; tempA ];
				b = [ b ; tempb ];
				Ae = [ Ae ; tempAe ];
				be = [ be ; tempbe ]; 

		    end

		    ibs.A = A;
		    ibs.b = b;
		    ibs.Ae = Ae;
		    ibs.be = be;
		    ibs.Dim = size(A,2);
		    
		    % If this is a closed loop set, then use the gain to close the loop!
		    if strcmp(ibs_settings.OpenLoopOrClosedLoop,'Closed')
		    	[tempA,tempb,tempAe,tempbe] = ibs.GetClosedLoopMatrices();
		    	%Save matrices
		    	ibs.A = tempA; ibs.b = tempb; ibs.Ae = tempAe; ibs.be = tempbe;
		    end

		    % Final Lang
		    FinalLang = KnowledgeSequence(end);

		    %%%%%%%%%%%%%%%%%%%%
			%% Create Outputs %%
			%%%%%%%%%%%%%%%%%%%%

		    ibs.AsPolyhedron =[];
		    ibs.t = t;
		    
		    ibs.System = lcsas;
		    ibs.KnowledgeSequence = KnowledgeSequence;

		end

		function [polyhedronOut] = ToPolyhedron(ibs)
			%Description:
			%	Converts InternalBehaviorSet to a Polyhedron() object from MPT3.

			polyhedronOut = Polyhedron( ...
				'A', ibs.A , 'b' , ibs.b , ...
				'Ae', ibs.Ae , 'be' , ibs.be ...
				);

			ibs.AsPolyhedron = polyhedronOut;
		end

		function [ polyhedronArrayOut ] = ToW_p(ibs)
			%Description:
			%	Projects the polyhedron onto the elements that represent 'w'.
			%

			%% Constants
			t = ibs.t;
			lcsas = ibs.System;
			L = lcsas.L;
			KnowledgeSequence = ibs.KnowledgeSequence;

			[ n_x , n_u , n_y , n_w , n_v ] = lcsas.Dimensions();

			%% Algorithm
			polyRepresentation = ibs.ToPolyhedron();
			polyRepresentation.computeVRep;

			polyhedronArrayOut = cell(L.cardinality(),1);

			LangAtEnd = KnowledgeSequence(end);
			for word_idx = 1:LangAtEnd.cardinality()
				projPolyhedron = polyRepresentation.projection(n_x*(t+1) + n_u*(t) + n_w*(t)*(word_idx-1) + [1:n_w*(t)]);
				projPolyhedron.computeHRep();
				
				[ ~ , wordIndexInL ] = L.contains(LangAtEnd.words{word_idx}); %Get index 
				polyhedronArrayOut{ wordIndexInL } = projPolyhedron;
				
			end

		end

		function [ ebs ] = ToExternalBehaviorSet(ibs)
			ebs = ExternalBehaviorSet(ibs.System,ibs.KnowledgeSequence);
		end

		function [ tf ] = CoversInputPolyhedron( ibs_array , poly_in )
			%CoversInputPolyhedron
			%Description:
			%	Checks whether or not the array of InternalBeliefSequence objects ibs_array
			%	spans enough of its space to cover the set poly_in.
			%	This is done using a simple optimization problem.
			%
			%Notes:
			%	In order to make MPT3 behave as expected, there is some bloating used. If the values in the problem are chosen randomly enough, this may lead to
			%	this function giving erroneous results.

			%% Input Processing

			if ~isa(poly_in,'Polyhedron')
				error(['Expected W_in to be a Polyhedron but instead it is of class ' class(W_in) '.' ])
			end

			%% Constants

			NumIBS = length(ibs_array);
			dim = poly_in.Dim;
		    
			if dim ~= ibs_array(1).Dim 
				error(['The dimension of the InternalBehaviorSet objects is ' num2str(ibs_array(1).Dim()) ' while the dimension of W is ' num2str(dim) '.' ])
			end

		    bloatFactor0 = 10^(-2);
            bloatPolyhedron = bloatFactor0 * Polyhedron('lb',-ones(1,dim),'ub',ones(1,dim));

		    verbosity = 0;

			%% Algorithm

			x = sdpvar(dim,1,'full');
			b = binvar(NumIBS,1);

			% Constrain x to be in poly_in
			x_in_p_constraint = [ poly_in.A * x <= poly_in.b ];
			if ~isempty(poly_in.Ae)
				x_in_p_constraint = x_in_p_constraint + [ poly_in.Ae * x == poly_in.be ];
			end

			% Constrain x so that if it is in any individual element of the polyunion
			% pu_in, then the binary flag is triggered.

			x_in_pu_constraint = [];
            ibs_array_as_polys = [];
			for set_index = 1:NumIBS
                ibs_array_as_polys = [ibs_array_as_polys;ibs_array(set_index).ToPolyhedron() + bloatPolyhedron ];
 				% temp_set = ibs_array_as_polys(set_index) * (1+bloatFactor0);

				x_in_pu_constraint = x_in_pu_constraint + [ implies( ibs_array_as_polys(set_index).A * x <= ibs_array_as_polys(set_index).b , b(set_index)) ];
				if ~isempty(ibs_array_as_polys(set_index).Ae)
					x_in_pu_constraint = x_in_pu_constraint + [ implies( ibs_array_as_polys(set_index).Ae * x == ibs_array_as_polys(set_index).be , b(set_index) ) ];
				end
			end

			% Create objective
			objective = sum(b);

			% Solve Optimization Flag
			ops = sdpsettings('verbose',verbosity,'debug',0);
			%ops = sdpsettings(ops,'solver','gurobi');
			
			optim0 = optimize(x_in_p_constraint+x_in_pu_constraint,objective,ops);

			if value(objective) == 0
				tf = false;
			else
				tf = true;
			end


		end

		function [ tf ] = CoversInputW( ibs_array , W_in , target_word )
			%CoversInputW
			%Description:
			%	Checks whether or not the array of InternalBeliefSequence objects ibs_array
			%	spans enough values of w to cover the set W_in.
			%	This is done using a simple optimization problem.
			%
			%Notes:
			%	In order to make MPT3 behave as expected, there is some bloating used. If the values in the problem are chosen randomly enough, this may lead to
			%	this function giving erroneous results.

			%% Input Processing

			if ~isa(W_in,'Polyhedron')
				error(['Expected W_in to be a Polyhedron but instead it is of class ' class(W_in) '.' ])
			end

			%% Constants

			NumIBS = length(ibs_array);
			W_dim = W_in.Dim;
			IBS_dim = ibs_array(1).Dim;

			System = ibs_array(1).System;
			t = ibs_array(1).t;
			[ n_x , n_u , n_y , n_w , n_v ] = System.Dimensions();

			if W_dim ~= (t*n_w) 
				error(['The dimension of the InternalBehaviorSet objects is ' num2str(IBS_dim) ' while the dimension of W is ' num2str(W_dim) '.' ])
			end

		    bloatFactor0 = 10^(-2);
		    verbosity = 1;

		    eps0 = 10^(-4);

		    cg = constr_gen(0);

			%% Algorithm
			x_cell_arr = {}; %x = sdpvar(ibs_array(1).Dim,1,'full');
			y_cell_arr = {};
			b = binvar(NumIBS,1);

			% Create Largest Possible Phi Set for This Word
			StarSequence = [ System.L ; repmat(Language(target_word),t,1) ];
			ibs_star = InternalBehaviorSet(System,StarSequence);

			phi = sdpvar(ibs_star.Dim,1,'full');
			phi_domain_constraint = [ ibs_star.A * phi <= ibs_star.b ] + [ ibs_star.Ae * phi == ibs_star.be ];

			% Constrain x so that if it is in any individual element of the polyunion
			% pu_in, then the binary flag is triggered.

			x_in_pu_constraint = []; w_matches_x_constraint = [];
            ibs_array_as_polys = [];
            objective_vec = [];
            word_selectors = {}; sf_arr = {};

            b_I = binvar(NumIBS,1);

            phi_implication_constraint = [];
			for set_index = 1:NumIBS
				
				%Get Current ibs, and make conditional constraint for phi if it is in this set.
				temp_ibs = ibs_array(set_index);

				if temp_ibs.KnowledgeSequence == StarSequence
					tf = true;
					return;
				end

				phi_implication_constraint = implies( [ temp_ibs.A * phi <= temp_ibs.b ] + [temp_ibs.Ae * phi == temp_ibs.be] , [b_I == 1] );

			end

			% Make sure that phi is in none of the given consistency sets
			not_included_in_any_constraint = [ sum(b_I) == 0 ];

			% Create objective
			objective = [];

			% Solve Optimization Flag
			ops = sdpsettings('verbose',verbosity,'debug',0);
			ops = sdpsettings(ops,'solver','gurobi');
			
			optim0 = optimize( ...
				phi_domain_constraint+phi_implication_constraint+not_included_in_any_constraint, ...
				objective, ...
				ops ...
			);

			if optim0.problem == 0
				%Problem is feasible. There exists a point that is not covered.
				tf = false;
			else
				%Problem is not feasible. There does not exist a point that is uncovered.
				disp(value(optim0.problem))
				tf = true;
			end


		end

		function [ constraints_out , implication_of_containment, word_selector , satisfaction_flag, x ] = CreateContainmentConstraintOnW(ibs,W_in)
			%Description:
			%	Creates a YALMIP constraint which determines that the vector w_in is contained/referenced
			%	by some element of ibs.
			%
			%Inputs:
			%	w - An sdpvar of appropriate dimension for the disturbance trajectory.
			%

			%% Input Checking

			% Make sure that only one internal behavior set was provided in the function call.
			if length(ibs) > 1
				error(['Call to CreateContainmentConstraintOnW() was done with an array of ' num2str(length(ibs)) ' InternalBehaviorSet objects. Expected only 1.' ])
			end

			%% Constants
			lcsas = ibs.System;
			t = ibs.t;
			[ n_x , n_u , n_y , n_w , n_v ] = lcsas.Dimensions();

			% if length(w_in) ~= n_w*t
			% 	error(['w input to CreateContainmentConstraintOnW() was expected to be of length ' num2str(n_w*t) ' but was instead of length ' num2str(length(w_in)) '.' ])
			% end

			KnowledgeSequence = ibs.KnowledgeSequence;
			LastLang = KnowledgeSequence(end);

			bloatFactor0 = 10^(-2);
			bloatedIBS = ibs.ToPolyhedron() + bloatFactor0 * Polyhedron('lb',-ones(1,ibs.Dim),'ub',ones(1,ibs.Dim));

			%% Algorithm

			word_selector = binvar(LastLang.cardinality(),1);
			satisfaction_flag = binvar(1,1);
			x = sdpvar(ibs.Dim,1,'full');

			constraints_out = [];

			%Create Containment Constraint
			% Add containment condition
			constraints_out = [ bloatedIBS.A * x <= bloatedIBS.b ];
			if ~isempty(ibs.Ae)
				constraints_out = constraints_out + [ bloatedIBS.Ae * x <= bloatedIBS.be ];
			end

			%Create implication which ties the conditions together

			selectionMatrices = ibs.R_T();

			for word_idx = 1:LastLang.cardinality()

				select_i = selectionMatrices{word_idx};

				%Add selection condition
				% constraints_out = constraints_out + [ implies( word_selector(word_idx) , w_in == select_i * x_in ) ];
				% constraints_out = constraints_out + [ implies( word_selector(word_idx) , w_in == select_i * x ) ];

				temp_containment_condition = [ W_in.A * select_i * x <= W_in.b ];
				if ~isempty(W_in.Ae)
					temp_containment_condition = temp_containment_condition + [ W_in.Ae * select_i * x <= W_in.be ];
				end

				constraints_out = constraints_out + [ implies( temp_containment_condition , satisfaction_flag ) ];

			end

			implication_of_containment = constraints_out;

			constraints_out = constraints_out + [ sum(word_selector) == 1 ]; %One of the word selector variables must be chosen.


		end

		function [ constraints_out , y_out , binvar_ibs_doesnt_contain ] = CreateNoncontainmentConditionOnW(ibs,w_in)
			%Description:
			%	For a given W_in, produces a constraint which says that the disturbance w_in
			%	does not appear in any tuple of the ibs phi = (x,u,w,x0).
			%
			%Usage:

			%% Constants
			selectionMatrices = ibs.R_T();
			cg = constr_gen(1);

			LastLang = ibs.KnowledgeSequence(end);

			%% Create Variables
			constraints_out = [];

			binvar_ibs_doesnt_contain = binvar(1,1);
			b = binvar(LastLang.cardinality(),1);

			bloatFactor0 = 10^(-2);
			bloatedIBS = ibs.ToPolyhedron() + bloatFactor0 * Polyhedron('lb',-ones(1,ibs.Dim),'ub',ones(1,ibs.Dim));

			%% Create Constraints

			H0 = bloatedIBS.A; h0 = bloatedIBS.b;
			if ~isempty(bloatedIBS.Ae)
				H0 = [H0;bloatedIBS.Ae;-bloatedIBS.Ae];
				h0 = [h0;bloatedIBS.be;-bloatedIBS.be];
			end

			H_f = {}; h_f = {};
			y_out = {};
			for word_index = 1:LastLang.cardinality()

				select_i = selectionMatrices{word_index};

				%Create a new constraint for full H
				H_f{word_index} = [H0;select_i;-select_i];
				h_f{word_index} = [h0;w_in;-w_in];

				[ temp_empty_constraint , temp_y ] = cg.polytope_is_empty_set( H_f{word_index} , h_f{word_index} );

				y_out{word_index} = temp_y;
				constraints_out = constraints_out + temp_empty_constraint; %[ implies(b(word_index),temp_nonempty_constraint) ];

				% y_out{word_index} = sdpvar(size(H_f{word_index},2),1);
				% constraints_out = constraints_out + [ iff( b(word_index) == 0 , H_f{word_index} * y_out{word_index} <= h_f{word_index} ) ]

			end

			%constraints_out = constraints_out + [ implies( binvar_ibs_doesnt_contain , [sum(b) == LastLang.cardinality() ] ) ];

		end

		function [ empty_flags ] = IsEmpty( ibs )
			%Description:
			%	This function identifies which of the InternalBehaviorSets are empty.
			%	If a given path's internal behavior set is empty, then we will mark that set appropriately.

			%% Input Checking

			if isscalar(ibs) && (~isnumeric([ibs.A,ibs.b]) || ~isnumeric([ibs.Ae,ibs.be]) )
				error('The matrices A,b, Ae, and be are expected to be numeric in order to check to see if this is an empty set.')
			end

			%% Variables

            empty_flags = logical.empty;
            
			%% Algorithm

			if isscalar(ibs)

				temp_poly = ibs.ToPolyhedron();
				empty_flags = temp_poly.isEmptySet;

			elseif isvector(ibs) && iscell(ibs)
                
				for ibs_index = 1:length(ibs)
					temp_single_ibs = ibs{ibs_index};
					empty_flags = [ empty_flags ; temp_single_ibs.IsEmpty() ];
				end

			elseif isvector(ibs) && ~iscell(ibs)

				for ibs_index = 1:length(ibs)
					temp_single_ibs = ibs(ibs_index);
					empty_flags = [ empty_flags ; temp_single_ibs.IsEmpty() ];
				end

			else
				error('Input to IsEmpty() must be a scalar or vector.')
			end

		end

		function [ selectionMatrices ] = R_T( ibs )
			%Description:
			%	Finds the matrix which selects the elements of a vector in the InternalBehaviorSet
			%	which correspond to the disturbances.
			%
			%Outputs
			%	selectionMatrices - A cell array containing n_w*t x ibs.Dim matrices.

			%% Constants
			lcsas = ibs.System;
			t = ibs.t;
			[ n_x , n_u , n_y , n_w , n_v ] = lcsas.Dimensions();
			FinalLang = ibs.KnowledgeSequence(end);

			%% Algorithm
			for word_idx = 1:FinalLang.cardinality()
				selectionMatrices{word_idx} = [ ...
					zeros(n_w*t, ibs.Dim - (FinalLang.cardinality() - word_idx + 1)*n_w*t - n_x*FinalLang.cardinality() ) , ...
					eye( n_w*t ) , ... 
					zeros(n_w*t, (FinalLang.cardinality() - word_idx) + n_x*FinalLang.cardinality() )  ...
					];
			end

		end

		function [ selectionMatrix ] = SelectExternalBehavior( ibs )
			%Description:
			%	Creates the matrix selectionMatrix that would retrieve the elements of a vector i in ibs such that 
			%	selectionMatrix * i  is the external behavior corresponding to i.

			%% Constants
			ibs_settings = ibs.ibs_settings;
			System = ibs.System;
			t = ibs.t;

			[ n_x , n_u , n_y , n_w , n_v ] = System.Dimensions();

			%% Algorithm 
			switch ibs_settings.fb_type
			case 'state'
				external_beh_dim = n_x*(t+1)+n_u*t;
			case 'output'
				external_beh_dim = n_y*(t+1)+n_u*t;
			otherwise
				error(['Unexpected fb_type given to internal_behavior_sets2containment_mat: ' fb_type ])
			end

			selectionMatrix = [ eye(external_beh_dim), zeros(external_beh_dim, ibs.Dim - external_beh_dim) ];

		end

		function [ selectionMatrices ] = SelectW( ibs )
			%Description:
			%	Finds a set of matrices which select each of the "w" matrix dimensions
			%	available in this Internal Behavior Set.
			%
			%Usage:
			%	selectionMatrices = ibs.SelectW();
			%

			%% Constants %%
			Dim = ibs.Dim;
			System = ibs.System;
			t = ibs.t;
			ibs_settings = ibs.ibs_settings;

			L = System.L;
			LastLang = ibs.KnowledgeSequence(end);

			[ n_x , n_u , n_y , n_w , n_v ] = System.Dimensions();

			%% Algorithm %%

			selectionMatrices = {};

			switch ibs_settings.fb_type
			case 'state'

				for ll_index = 1:LastLang.cardinality()
					temp_word = LastLang.words{ll_index};
					[ ~ , tw_index ] = L.contains(temp_word);

					selectionMatrices{ll_index} = [ ...
						zeros(n_w*t,n_x*(t+1) + n_u*t + (tw_index-1)*n_w*t ) , ...
						eye(n_w*t) , ...
						zeros( n_w*t , Dim-n_x*(t+1)-n_u*t-tw_index*n_w*t) ...
						];

				end
			case 'output'

				error(['The output feedback version of SelectW() is not working yet!'])

			otherwise
				error(['Unexpected Feedback Type: ' ibs_settings.fb_type ])
			end


		end

		function [ selectionMatrices ] = SelectX0( ibs )
			%Description:
			%	Finds a set of matrices which select each of the "x0" components of an element
			%	from the internal behavior set ibs.
			%
			%Usage:
			%	selectionMatrices = ibs.SelectX0()

			%% Constants %%
			Dim = ibs.Dim;
			System = ibs.System;
			t = ibs.t;
			ibs_settings = ibs.ibs_settings;

			FirstLang = ibs.KnowledgeSequence(1);
			LastLang = ibs.KnowledgeSequence(end);

			[ n_x , n_u , n_y , n_w , n_v ] = System.Dimensions();

			%% Algorithm %%

			selectionMatrices = {};

			switch ibs_settings.fb_type
			case 'state'
				% Matrices will all be the same because of how the set was defined.

				for word_index = 1:LastLang.cardinality()
					selectionMatrices{word_index} = [ ...
						zeros(n_x,n_x*(t+1) + n_u*t + FirstLang.cardinality()*n_w*t ) , ...
						eye(n_x) ...
						];

				end
			case 'output'

				error(['The output feedback version of SelectW() is not working yet!'])

			otherwise
				error(['Unexpected Feedback Type: ' ibs_settings.fb_type ])
			end
		end

	end

end