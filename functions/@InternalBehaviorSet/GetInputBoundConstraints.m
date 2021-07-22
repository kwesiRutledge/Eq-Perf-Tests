function [ constraints , dual_vars ] = GetInputBoundConstraints( ibs )
	%Description:
	%	Creates a containment constraint based on the internal behavior set ibs and the
	%	desired bound on the input set.
	%	The time horizon for reachability is assumed to be given by the system's language.
	%
	%Usage:
	%	[ constraints , dual_vars ] = ibs.GetInputBoundConstraints()

	%% Input Processing %%

	if ~strcmp(ibs.ibs_settings.fb_type,'state')
		error(['GetReachabilityConstraints() was designed for state feedback only!'])
	end

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
	U_T = 1;
	for tau = 1:T
		U_T = U_T * U;
	end

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
	SelectWMatrices = ibs.SelectW();
	selectW1 = SelectWMatrices{1};

	selectWAll = [ 	selectW1, zeros(n_w*(T-1),n_w) ; 
					[ zeros(n_w,ibs.Dim) , eye(n_w) ] ];

	% Create Contraints

	% [ dual_vars , constraints ] = cg.create_sadraddini_AH_inclusion_constr( ...
	% 	zeros(n_w*T,1) , selectWAll   , H , h , ...
	% 	zeros(n_w*T,1) , [eye(n_w*T)] , U_T.A*ibs.K , U_T.b - U_T.A*ibs.k );

	[ dual_vars , constraints ] = cg.get_H_polyt_inclusion_constr( H , h , U_T.A*ibs.K*selectWAll , U_T.b - U_T.A*ibs.k );

end