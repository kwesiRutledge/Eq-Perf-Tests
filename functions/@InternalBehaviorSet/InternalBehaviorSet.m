classdef InternalBehaviorSet
	%Description:
	%	The internal behavior set associated with a given LCSAS.
	%
	%To-do:
	%	- Add second constructor which does what get_closed_loop_consistent_internal_behavior_set_matrices.m is doing now.

	properties
		System;
		Polyhedron;
		t;
		KnowledgeSequence;

		% H Representation
		A;
		b;
		Ae;
		be;

		% Other Things
		Dim;
	end

	methods
		function [ibs] = InternalBehaviorSet(varargin)
			%internal_behavior_set.m
			%Description:
			%	Finds a polyhedron that describes what pairs of state, input, and disturbance sequences are compatible/feasible from ALL
			%	switching sequences defined by the sequence KnowledgeSequence.
			%
			%	The KnowledgeSequence is a vector describing the following sequence of estimates about the true state of the system;
			%		KnowledgeSequence = [ L_{0} , L_{1} , ... , L_{T-1} ]
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
			%					{ [y]  | }
			%					{ [u]  | }
			%		C(\sigma) = { [w]  | }
			%					{ [v]  | }
			%					{ [x0] | }
			%					{ [x]  | }
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
			[ lcsas , KnowledgeSequence , ibs_settings ] = ip_InternalBehaviorSet( varargin{:} );

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
			T = length(lcsas.L.words{1});
		    
		    P_u = lcsas.U;

		    [ n_x , n_u , n_y , n_w , n_v ] = lcsas.Dimensions();

			%% Algorithm
		        
			% Debugging Strategy
			if ibs_settings.debug > 0
				disp('Debugging info...')
			end

		    % Create Internal Behavior Set Using Each of the ClairvoyantSequences
		        
		    ib_components_at_ = {};
		    % ib_components_at_{1} = lcsas0.X0;

		    for t = 1:KS_len %KnowledgeSequence estimates L_0, ... , L_{T-2}
		        L_t = KnowledgeSequence(t); 

		        %Collect the ib_at_t sets for each word in L_t
		        ib_components_at_{t} = {};
		        for word_idx = 1:L_t.cardinality()
		            L_t_projected = Language(L_t.words{word_idx});
		            temp_ib_sets(word_idx) = internal_behavior_set_at_t( lcsas , t , L_t_projected , ibs_settings );

		            % Transform the set into the dimension described by max_card.
		            if word_idx == 1
		                ib_components_at_{t} = modify_ibs_with_respect_to_T( lcsas , KnowledgeSequence , temp_ib_sets(word_idx) , L_t_projected , t , KS_len );
		            else
		                ib_components_at_{t} = ib_components_at_{t}.intersect( ...
		                    modify_ibs_with_respect_to_T( lcsas , KnowledgeSequence , temp_ib_sets(word_idx) , L_t_projected , t , KS_len ) ...
		                );	
		            end
		        end
		    end

		    ib_at_T = ib_components_at_{1};
		    for t = 2:length(ib_components_at_)
		        ib_at_T = ib_at_T.intersect(ib_components_at_{t});
		    end
		    
		    % Final Lang
		    FinalLang = KnowledgeSequence(end);

		    %%%%%%%%%%%%%%%%%%%%
			%% Create Outputs %%
			%%%%%%%%%%%%%%%%%%%%

		    ibs.Polyhedron = ib_at_T ;
		    ibs.t = KS_len;
		    ibs.Dim = n_x * (ibs.t+1) + n_u*ibs.t + n_w*ibs.t*FinalLang.cardinality() + n_x*FinalLang.cardinality();
		    
		    ibs.A = ib_at_T.A;
		    ibs.b = ib_at_T.b;
		    ibs.Ae = ib_at_T.Ae;
		    ibs.be = ib_at_T.be;
		    
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
		end

		function [ polyhedronArrayOut ] = ToW(ibs)
			%Description:
			%	Projects the polyhedron onto the elements that represent 'w'.
			%

			%% Constants
			t = ibs.t;
			lcsas = ibs.System;
			KnowledgeSequence = ibs.KnowledgeSequence;

			[ n_x , n_u , n_y , n_w , n_v ] = lcsas.Dimensions();

			%% Algorithm
			polyRepresentation = ibs.ToPolyhedron();

			polyhedronArrayOut = [];
			LangAtEnd = KnowledgeSequence(end);
			for word_idx = 1:LangAtEnd.cardinality()
				projPolyhedron = polyRepresentation.projection(n_x*(t+1) + n_u*(t) + n_w*(t)*(word_idx-1) + [1:n_w*(t)]);
				polyhedronArrayOut = [ polyhedronArrayOut ; projPolyhedron ];
			end

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

		function [ tf ] = CoversInputW( ibs_array , W_in )
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

			lcsas = ibs_array(1).System;
			t = ibs_array(1).t;
			[ n_x , n_u , n_y , n_w , n_v ] = lcsas.Dimensions();

			if W_dim ~= (t*n_w) 
				error(['The dimension of the InternalBehaviorSet objects is ' num2str(t*n_w) ' while the dimension of W is ' num2str(W_dim) '.' ])
			end

		    bloatFactor0 = 10^(-2);
		    verbosity = 1;

		    eps0 = 10^(-4);

		    cg = constr_gen();

			%% Algorithm

			w = sdpvar(W_dim,1,'full');
			x_cell_arr = {}; %x = sdpvar(ibs_array(1).Dim,1,'full');
			y_cell_arr = {};
			b = binvar(NumIBS,1);

			% Constrain x to be in W_in
			x_in_p_constraint = [ W_in.A * w <= W_in.b ];
			if ~isempty(W_in.Ae)
				x_in_p_constraint = x_in_p_constraint + [ W_in.Ae * w == W_in.be ];
			end

			% Constrain x so that if it is in any individual element of the polyunion
			% pu_in, then the binary flag is triggered.

			x_in_pu_constraint = []; w_matches_x_constraint = [];
            ibs_array_as_polys = [];
            objective_vec = [];
            word_selectors = {}; sf_arr = {};
			for set_index = 1:NumIBS
				
				% Create Constraint Component
				% [ temp_containment_condition , temp_implication , temp_word_selector , satisfaction_flag , x_out ] = ibs_array(set_index).CreateContainmentConstraintOnW( W_in );
				% word_selectors{end+1} = temp_word_selector;
				% x_cell_arr{end+1} = x_out;
				% sf_arr{end+1} = satisfaction_flag;
				% x_in_pu_constraint = x_in_pu_constraint + temp_containment_condition; %temp_implication + [ sum(word_selectors{end}) == b(set_index) ];

				% temp_ibs = ibs_array(set_index);
				% temp_H = temp_ibs.A; temp_h = temp_ibs.b;
				% if ~isempty(temp_ibs.Ae)
				% 	temp_H = [temp_H;temp_ibs.Ae;-temp_ibs.Ae];
				% 	temp_h = [temp_h;temp_ibs.be;-temp_ibs.be];
				% end

				temp_ibs = ibs_array(1);
				[ temp_empty_constraint , temp_y , doesnt_contain ] = temp_ibs.CreateNoncontainmentConditionOnW( w );
				y_cell_arr{end+1} = temp_y;

				x_in_pu_constraint = x_in_pu_constraint + temp_empty_constraint; %+ [ b(set_index) == doesnt_contain ];

				%Add to objective;
				%objective_vec = [objective_vec;sf_arr{end}];

			end

			% Create objective
			%objective = sum(b);
			%objective = min(objective_vec);
			% objective = -sum(b);
			objective = [];

			% Solve Optimization Flag
			ops = sdpsettings('verbose',verbosity,'debug',0);
			ops = sdpsettings(ops,'solver','gurobi');
			
			optim0 = optimize(x_in_p_constraint+x_in_pu_constraint,objective,ops);

			% if value(objective) == NumIBS
			% 	% disp(value(b))
			% 	% disp(value(w))
			% 	disp(value(objective_vec))

			% 	tf = false;
			% else
			% 	tf = true;
			% end

			if optim0.problem == 0
				%Problem is feasible. There exists a point that is not covered.
				disp(value(w))
				disp(value(w(2)))

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

	end

end

%% Class-Related Functions

function [ lcsas0 , KnowledgeSequence , ibs_settings ] = ip_InternalBehaviorSet( varargin )
	%Description:
	%	Input processing function for internal_behavior_set().
	%	Verifies that there:
	%	- are at least four inputs,
	%	- LCSAS input checks,

	%% Check Argument Numbers
	expected_number_of_arguments = 2;

	if nargin < expected_number_of_arguments
		error(['Expect at least four arguments.'])
	end

	%% First Required Input: lcsas
	lcsas0 = varargin{1};
	lcsas0.check('X0','U'); %Verify that X0 and U exist and are both Polyhedron objects.

	%% Second Required Input: KnowledgeSequence
	KnowledgeSequence = varargin{2};
	if ~isa(KnowledgeSequence,'Language')
		error(['KnowledgeSequence  is not of type ''Language''. Instead it is of type ' class(KnowledgeSequence) '.'])
	end

	if ~isvector(KnowledgeSequence)
		error(['KnowledgeSequence must be a vector or scalar of Language objects.'])
	end

	%% Define Defaults
	ibs_settings = struct('fb_type','state','reduce_flag',false,'debug',0,'ReturnEarly',false, ...
							'A', [] , 'b' , [] , 'Ae' , [] , 'be' , [] , ...
							'Dim', -1 );

	all_fb_types = {'state','output'};

	%% Check the additional inputs
	argument_index = expected_number_of_arguments+1;
	while argument_index <= nargin
		switch varargin{argument_index}
			case 'fb_type'
				ibs_settings.fb_type = varargin{argument_index+1};
				if ~any(strcmp(ibs_settings.fb_type,all_fb_types))
					error(['Unexpected fb_type value: ' ibs_settings.fb_type ])
				end
				argument_index = argument_index + 2;
			case 'reduce_flag'
				ibs_settings.reduce_flag = varargin{argument_index+1};
				argument_index = argument_index + 2;
			case 'ReturnEarly'
				ibs_settings.ReturnEarly = varargin{argument_index+1};
				argument_index = argument_index + 2;
			case 'A'
				ibs_settings.A = varargin{argument_index+1};
				argument_index = argument_index + 2;
			case 'b'
				ibs_settings.b = varargin{argument_index+1};
				argument_index = argument_index + 2;
			case 'Ae'
				ibs_settings.Ae = varargin{argument_index+1};
				argument_index = argument_index + 2;
			case 'be'
				ibs_settings.be = varargin{argument_index+1};
				argument_index = argument_index + 2;
			otherwise
				error(['Unexpected input to internal_behavior_set: ' varargin{argument_index} ])
		end
	end

	% Check inputs when the ReturnEarly flag is raised.
	if ibs_settings.ReturnEarly
		if isempty(ibs_settings.A) && isempty(ibs_settings.Ae)
			error(['The InternalBehaviorSet constructed using ReturnEarly must contain a form of H constraint.'])
		end

		tempDimArray = [ size(ibs_settings.A,2) ; size(ibs_settings.Ae,2) ];
		if all(tempDimArray > 0) && ( tempDimArray(1) ~= tempDimArray(2) )
			error(['The dimension of the set according to A is ' num2str(size(ibs_settings.A,2)) ' while the dimension of the set according to Ae is ' num2str(size(ibs_settings.Ae,2)) '. Please adjust your input matrices.'  ])
		end 

		ibs_settings.Dim = max(tempDimArray);

	end
end

function [ ib_at_t ] = internal_behavior_set_at_t( lcsas0 , t , L_t , ibs_settings )
	%Description:
	%	Creates the set of behaviors that are consistent with the hypothesis L_t at time t.
	%	This set does not consider hypotheses at previous times and is thus a sometimes over approximation
	%	of what can be expected.
	%
	%Notes:
	%	It's possible that this can be simplified. We are enforcing a lot of constraints that don't need to be enforced.
	%	i.e. the set is enforcing constraints also on times t-1,t-2, etc. but it does not need to.
	%
	%	Each trajectory contained in the internal behavior set starts from t=0. So, the internal behavior_set_at_t=1
	%	will contain 
	%
	%			[ x0 ]
	%		x = [ x1 ] , u = [ u0 ] , w = [ w0 ]
	%

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	[ n_x , n_u , n_y , n_w , n_v ] = lcsas0.Dimensions();

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	% Create large disturbance set from Cartesian products
	P_wT = 1; P_x0_L = 1;
	for word_idx = 1:L_t.cardinality()
		for symb_idx = 1:t
			P_wT = P_wT * lcsas0.Dyn( L_t.words{word_idx}(symb_idx) ).P_w;
		end
		%Initial Condition Set for Each Word in L
		P_x0_L = P_x0_L * lcsas0.X0;
    end
    
    % Created Disturbance Sets
    P_uT = 1;
    for t_idx = 1:t
        P_uT = P_uT * lcsas0.U;
    end

    % Create mpc matrices for each word in the language L
	Hc = {}; Sc = {}; Jc = {}; fc = {}; Cc = {}; Bwc = {}; Cvc = {};
	for word_ind = 1:length(L_t.words)
		[Hc{word_ind},Sc{word_ind},Cc{word_ind},Jc{word_ind},fc{word_ind},Bwc{word_ind},Cvc{word_ind}] = lcsas0.get_mpc_matrices('word',L_t.words{word_ind}(1:t));
	end

	H_block = []; S_block = []; J_block = []; f_block = [];
	I_blockx = []; I_blockx2 = [];  I_blocky = [];
	C_block = []; Cv_block = [];
	for word_ind = 1:length(L_t.words)
		H_block(end+[1:size(Hc{word_ind},1)],end+[1:size(Bwc{word_ind},2)]) = Hc{word_ind}*Bwc{word_ind};
		S_block(end+[1:size(Sc{word_ind},1)],[1:size(Sc{word_ind},2)]) = Sc{word_ind};
		J_block(end+[1:size(Jc{word_ind},1)],end+[1:size(Jc{word_ind},2)]) = Jc{word_ind};
		f_block(end+[1:size(Hc{word_ind}*fc{word_ind},1)],1) = Hc{word_ind}*fc{word_ind};

		I_blockx(end+[1:n_x*(t+1)],[1:n_x*(t+1)]) = eye(n_x*(t+1));
		I_blocky(end+[1:n_y*(t+1)],[1:n_y*(t+1)]) = eye(n_y*(t+1));

	end

	%% Constructing the Sets
	if strcmp(ibs_settings.fb_type,'state')
	    
		P_eta = P_uT * P_wT * P_x0_L;

		%Create the set of feasible (x,u,w,x0) tuples
		ib_at_t = Polyhedron(	'A',[zeros(size(P_eta.A,1),n_x*(t+1)),P_eta.A],'b',P_eta.b, ...
								'Ae',[-I_blockx, S_block, H_block, J_block],'be',-f_block );

		% Add Equality Constraints Between Words in Internal Behavior Sets
		for word_ind = 2:L_t.cardinality()
			Ltc = L_t.cardinality();
			matching_x_block = [ eye(n_x*(t+1)) , zeros( n_x*(t+1) , (-2+word_ind)*n_x*(t+1) ) , -eye(n_x*(t+1)) , zeros(n_x*(t+1) , (Ltc-word_ind)*n_x*(t+1) ) ];
			ib_at_t.Ae = [ 	ib_at_t.Ae ;
							matching_x_block , zeros(n_x*(t+1),size(S_block,2)) , zeros(n_x*(t+1),size(H_block,2)) , zeros(n_x*(t+1),size(J_block))  ];

			ib_at_t.be = [ ib_at_t.be ; zeros(n_x*(t+1),1) ];
		end

	else strcmp(ibs_settings.fb_type,'output')

		%Also introduce the measurement disturbance into the equation
		P_vT = 1; 
		for word_idx = 1:length(L.words)
			for symb_idx = 1:t+1
				P_vT = P_vT * lcsas.Dyn( L.words{word_idx}(symb_idx) ).P_v;
			end
    	end

    	% Introduce Sets
    	for word_ind = 1:L.cardinality()
    		C_block(end+[1:n_y*(t+1)],end+[1:n_x*(t+1)]) = [Cc{word_ind} ; zeros(n_y,n_x*t), lcsas.Dyn( L.words{word_ind}(t+1) ).C ];
			Cv_block(end+[1:n_y*(t+1)],end+[1:n_v*(t+1)]) = [Cvc{word_ind},zeros(size(Cvc{word_ind},1),n_v);zeros(n_y,size(Cvc{word_ind},2)), lcsas.Dyn( L.words{word_ind}(t+1) ).C_v ];
			I_blockx2(end+[1:n_x*(t+1)],end+[1:n_x*(t+1)]) = eye(n_x*(t+1));
		end

    	P_eta = P_uT * P_wT * P_vT * P_x0_L;

    	%Create the set of feasible (x,u,w,x0) tuples
    	ib_at_t = Polyhedron(	'A',[zeros(size(P_eta.A,1),n_y*(t+1)),P_eta.A,zeros(size(P_eta.A,1),length(L.words)*n_x*(t+1))],'b',P_eta.b, ...
    							'Ae',[zeros(size(S_block,1),size(I_blocky,2)),S_block, H_block, zeros(size(S_block,1),size(Cv_block,2)), J_block, -I_blockx2; ...
    								  I_blocky, zeros(size(I_blocky,1),size(S_block,2)+size(H_block,2)), -Cv_block , zeros(size(I_blocky,1),size(J_block,2)) , -C_block ], ...
    							'be', [-f_block;zeros(size(I_blocky,1),1)] );
    end

    if (~ib_at_t.isEmptySet) && ibs_settings.reduce_flag
		ib_at_t.minHRep; %Make sure to do this to simplify some of the future projections.
	end

end

function [ ibs_extended ] = modify_ibs_with_respect_to_T( lcsas0 , KnowledgeSequence , ibs_in , L_in , t , T )
	%Description:
	%	Modifies the internal behavior set defined for a time t and map it into the dimension
	%	of the expected behaviors at time T.
	%

	%% Input Checking

	% Verify that the input language has a single word in it.
	if L_in.cardinality() ~= 1
		error(['modify_ibs_with_respect_to_T() expects a language with only one word inside of it.'])
	end

	%% Constants

	[ n_x , n_u , n_y , n_w , n_v ] = lcsas0.Dimensions();

    % trimmedKnowledgeSequence = KnowledgeSequence([2:end]);
    
	[ max_card , max_card_index ] = KnowledgeSequence.find_maximum_cardinality_in_sequence();

	ibs_e_dim = n_x * (T+1) + n_u*T + n_w*T*max_card + n_x*max_card;

	%% Algorithm

	%Find the location of the input Language L_in relative to the language at max_card_t
	L_mct = KnowledgeSequence(max_card_index);
	[tf,L_in_index] = L_mct.contains(L_in.words{1});
 
	if ~tf
		error(['The input Language L_in for modify_ibs_with_respect_to_T() is not contained within L_mct.'])
	end

	L_in_index_as_binary = zeros(1,L_mct.cardinality());
	L_in_index_as_binary( L_in_index ) = 1;

	A_Prefactor = [ ...
		eye( n_x * (t+1) ), zeros( n_x*(t+1) , ibs_e_dim - n_x * (t+1) ) ;
		zeros( n_u*(t) , n_x*(T+1) ), eye(n_u*t) , zeros( n_u*(t) , ibs_e_dim - n_x*(T+1) - n_u*t ) ;
		zeros( n_w*(t) , n_x*(T+1) + n_u*T ) , kron( L_in_index_as_binary , [ eye(n_w*t) , zeros(n_w*t,n_w*(T-t)) ] ) , zeros(n_w*t, ibs_e_dim - n_x*(T+1) - n_u*T - n_w*T*max_card) ;
		zeros( n_x , ibs_e_dim - n_x*max_card), kron( L_in_index_as_binary , eye(n_x) ) ...
		];

	ibs_extended = Polyhedron( ...
		'A' , ibs_in.A *A_Prefactor, 'b' , ibs_in.b , ...
		'Ae', ibs_in.Ae*A_Prefactor, 'be', ibs_in.be );

end 