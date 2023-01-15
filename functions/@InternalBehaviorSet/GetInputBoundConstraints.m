function [ constraints , dual_vars ] = GetInputBoundConstraints( varargin )
	%Description:
	%	Creates a containment constraint based on the internal behavior set ibs and the
	%	desired bound on the input set.
	%	The time horizon for reachability is assumed to be given by the system's language.
	%
	%Usage:
	%	[ constraints , dual_vars ] = ibs.GetInputBoundConstraints()
	%	[ constraints , dual_vars ] = ibs.GetInputBoundConstraints('Relaxation',true)

	%% Input Processing %%
	[ ibs , use_cl_or_ol_A ] = ip_GetInputBoundConstraints(varargin{:});

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	A = ibs.A; b = ibs.b; Ae = ibs.Ae; be = ibs.be; %Get Polytope Matrices

	System = ibs.System;
	L = System.L;

	T = length(L.words{1});
	[ n_x , n_u , n_y , n_w , n_v ] = System.Dimensions();

	KnowledgeSequence = ibs.KnowledgeSequence;
	LastLang = KnowledgeSequence(end);

	cg = constr_gen(0);

	U = System.U;

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	switch use_cl_or_ol_A
	case 'A_ol'
		%Use Open Loop Relaxation
		[constraints , dual_vars] = GetRelaxedInputBoundConstraints(ibs);
	case 'A_cl'
		%Use Full Closed Loop Matrices to Create Constraint
		[constraints , dual_vars] = GetCompleteInputBoundConstraints(ibs);
	otherwise
		error(['Unexpected value for use_cl_or_ol_A in GetInputBoundConstraints().'])
	end

end

function [ ibs , use_cl_or_ol_A ] = ip_GetInputBoundConstraints(varargin)
	%Description:
	%	Creates the potential inputs from the call to ip.

	%% Get IBS %%

	ibs = varargin{1};

	if ~strcmp(ibs.ibs_settings.fb_type,'state')
		error(['GetReachabilityConstraints() was designed for state feedback only!'])
	end

	%% Check To See If Any Desired Settings Were Given %%
	use_cl_or_ol_A = 'A_ol';

	argument_index = 2;
	while argument_index <= nargin
		switch varargin{argument_index}
		case 'Use A_cl or A_ol?'
			use_cl_or_ol_A = varargin{argument_index+1};
			argument_index = argument_index + 2;
		otherwise
			error(['Unexpected input to GetInputBoundConstraints(): ' varargin{argument_index} ])
		end
	end


end 

function [ constraints , dual_vars ] = GetRelaxedInputBoundConstraints(ibs)
	%Description:
	%	Removes the bilinearity in the containment condition by creating the constraints
	%	with respect to the OPEN LOOP version of ibs.

	%% Constants

	ibs_ol = InternalBehaviorSet(ibs.System,ibs.KnowledgeSequence); %Get Open Loop Version of ibs

	A = ibs_ol.A; b = ibs_ol.b; 
	Ae = ibs_ol.Ae; be = ibs_ol.be; %Get Polytope Matrices

	System = ibs.System;
	L = System.L;

	T = length(L.words{1});
	[ n_x , n_u , n_y , n_w , n_v ] = System.Dimensions();

	KnowledgeSequence = ibs.KnowledgeSequence;
	LastLang = KnowledgeSequence(end);

	cg = constr_gen(0);

	U = System.U;

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	% Create T repeated version of U
	U_Tm1 = 1;
	for tau = 1:T-1
		U_Tm1 = U_Tm1 * U;
	end
	U_T = U_Tm1 * U;

	constraints = [];
	dual_vars = {};
	for word_index = 1: LastLang.cardinality()

		%Use first word in Last Language (of Knowledge Sequence) to Help Create H
		word1 = LastLang.words{1};
		lastSymbol1 = word1(end);

		[~,word1_index] = L.contains(word1);

		W_final = System.Dyn(lastSymbol1).P_w;

		% Create the Set Which Will Represent The Possible Sequence of Disturbances

		H = [ A ; Ae ; -Ae];
		H = [ H , zeros(size(H,1),n_w) ; zeros(size(W_final.A,1),ibs.Dim) , W_final.A ];

		h = [ b ; be ; -be ; W_final.b];

		if ~isempty(W_final.Ae)
			H = [ H ; zeros(size(W_final.Ae,1),ibs.Dim) , W_final.Ae ; zeros(size(W_final.Ae,1),ibs.Dim) , -W_final.Ae];
			h = [ h ; W_final.be ; -W_final.be ];
		end

		% Create Input Set Matrices

		% Create Containment Based on Full Length Disturbances from H
		selectW1 = ibs.SelectW();

		selectWAll = [ 	selectW1, zeros(n_w*(T-1),n_w) ; 
						[ zeros(n_w,ibs.Dim) , eye(n_w) ] ];

		% Create Contraints

		[ dual_vars , constraints ] = cg.get_H_polyt_inclusion_constr( H , h , U_T.A*ibs.K*selectWAll , U_T.b - U_T.A*ibs.k );

	end

end

function [ constraints , dual_vars ] = GetCompleteInputBoundConstraints(ibs)
	%Description:
	%	

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	A = ibs.A; b = ibs.b; Ae = ibs.Ae; be = ibs.be; %Get Polytope Matrices

	System = ibs.System;
	L = System.L;

	T = length(L.words{1});
	[ n_x , n_u , n_y , n_w , n_v ] = System.Dimensions();

	KnowledgeSequence = ibs.KnowledgeSequence;
	LastLang = KnowledgeSequence(end);

	cg = constr_gen(0);

	U = System.U;

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	% Create T repeated version of U
	U_Tm1 = 1;
	for tau = 1:T-1
		U_Tm1 = U_Tm1 * U;
	end
	U_T = U_Tm1 * U;

	%Use first word in Last Language (of Knowledge Sequence) to Help Create H
	constraints = [];
	dual_vars = {};
	for word_index = 1: LastLang.cardinality()
		word1 = LastLang.words{word_index};
		lastSymbol1 = word1(end);

		[ ~ , word1_index ] = L.contains(word1);
		W_final = System.Dyn(lastSymbol1).P_w;

		% Create the Set Which Will Represent The Possible Sequence of Disturbances

		H = [ A ; Ae ; -Ae];
		H = [ H , zeros(size(H,1),n_w) ; zeros(size(W_final.A,1),ibs.Dim) , W_final.A ];

		h = [ b ; be ; -be ; W_final.b];

		if ~isempty(W_final.Ae)
			H = [ H ; zeros(size(W_final.Ae,1),ibs.Dim) , W_final.Ae ; zeros(size(W_final.Ae,1),ibs.Dim) , -W_final.Ae];
			h = [ h ; W_final.be ; -W_final.be ];
		end

		% Create Input Set Matrices

		% Create Containment Based on Full Length Disturbances from H
		selectW1 = ibs.SelectW();

		selectWAll = [ 	selectW1, zeros(n_w*(T-1),n_w) ; 
						[ zeros(n_w,ibs.Dim) , eye(n_w) ] ];

		% Create Contraints
		[ dual_vars{1+2*(word_index-1)} , temp_constraints ] = cg.get_H_polyt_inclusion_constr( ...
			H , h , U_T.A * ibs.K * selectWAll , U_T.b - U_T.A * ibs.k ...
		);
		constraints = constraints + temp_constraints;

		% [ dual_vars{2+2*(word_index-1)} , temp_constraints ] = cg.get_H_polyt_inclusion_constr( H , h , U.A*ibs.K(end-[1:n_u],:)*selectWAll , U.b - U.A*ibs.k(end-[1:n_u]) );
		% constraints = constraints + temp_constraints;
	end

end