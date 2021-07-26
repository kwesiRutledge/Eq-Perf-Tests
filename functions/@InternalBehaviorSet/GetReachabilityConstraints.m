function [ constraints , dual_vars ] = GetReachabilityConstraints( varargin )
	%Description:
	%	Creates a reachability constraint based on the internal behavior set ibs and the
	%	target set.
	%	The time horizon for reachability is assumed to be given by the system's language.
	%
	%Usage:
	%	[ constraints , dual_vars ] = ibs.GetReachabilityConstraints( X_Target )
	%	[ constraints , dual_vars ] = ibs.GetReachabilityConstraints( X_Target , 'Relaxation' , true )

	%% Input Processing %%

	[ ibs , X_Target , relax ] = ip_GetReachabilityConstraints(varargin{:});

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	System = ibs.System;
	L = System.L;

	T = length(L.words{1});
	[ n_x , n_u , n_y , n_w , n_v ] = System.Dimensions();

	KnowledgeSequence = ibs.KnowledgeSequence;
	LastLang = KnowledgeSequence(end);

	cg = constr_gen(0);

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	%Use first word in Last Language (of Knowledge Sequence) to Help Create H
	word1 = LastLang.words{1};
	lastSymbol1 = word1(end);

	W_final = System.Dyn(lastSymbol1).P_w;

	% Create the Set Which Will Represent The Possible Sequence of Disturbances
	if relax
		ibs_ol = InternalBehaviorSet(ibs.System,ibs.KnowledgeSequence); %Get Open Loop Version of ibs

		A = ibs_ol.A; b = ibs_ol.b; 
		Ae = ibs_ol.Ae; be = ibs_ol.be; %Get Polytope Matrices

	else
		A = ibs.A; b = ibs.b; Ae = ibs.Ae; be = ibs.be; %Get Polytope Matrices
	end

	H = [ A ; Ae ; -Ae];
	H = [ H , zeros(size(H,1),n_w) ; zeros(size(W_final.A,1),ibs.Dim) , W_final.A ];

	h = [ b ; be ; -be ; W_final.b];

	if ~isempty(W_final.Ae)
		H = [ H ; zeros(size(W_final.Ae,1),ibs.Dim) , W_final.Ae ; zeros(size(W_final.Ae,1),ibs.Dim) , -W_final.Ae];
		h = [ h ; W_final.be ; -W_final.be ];
	end

	% Create Containment Based on Full Length Disturbances from H
	SelectWMatrices = ibs.SelectW();
	SelectX0Matrices = ibs.SelectX0();

	selectW1 = SelectWMatrices{1};
	selectX0 = SelectX0Matrices{1};

	selectWAndX0 = [ 	selectW1, zeros(n_w*(T-1),n_w) ; 
						[ zeros(n_w,ibs.Dim) , eye(n_w) ]  ;  
						selectX0 , zeros(n_x,n_w) ];

	% Create Target Matrices
	[S_w,S_u,S_C,S_x0,S_k,S_Bw,S_Cv] = System.get_mpc_matrices('word',word1); %Get MPC Matrices
	G = [ 	S_w*S_Bw+S_u*ibs.K , ...
			S_x0 ];

	selectFinalState = [ zeros(n_x,n_x*T) , eye(n_x) ];

	H_T = X_Target.A * selectFinalState * G * selectWAndX0;
	h_T = X_Target.b - X_Target.A * selectFinalState * ( S_w * S_k + S_u * ibs.k);

	% [ dual_vars , constraints ] = cg.create_sadraddini_AH_inclusion_constr( ...
	% 	zeros(n_w*T+n_x,1) , selectWAndX0   , H , h , ...
	% 	zeros(n_w*T+n_x,1) , [eye(n_w*T+n_x)] , H_T , h_T );

	[ dual_vars , constraints ] = cg.get_H_polyt_inclusion_constr( H , h , H_T, h_T );


end

function [ ibs , X_Target , relax ] = ip_GetReachabilityConstraints(varargin)
	%Description:
	%	Creates the potential inputs from the call to ip.

	%% Get IBS %%

	ibs = varargin{1};
	X_Target = varargin{2};

	if ~strcmp(ibs.ibs_settings.fb_type,'state')
		error(['GetReachabilityConstraints() was designed for state feedback only!'])
	end

	%% Check To See If Any Desired Settings Were Given %%
	relax = false;

	argument_index = 3;
	while argument_index <= nargin
		switch varargin{argument_index}
		case 'Relaxation'
			relax = varargin{argument_index+1};
			argument_index = argument_index + 2;
		otherwise
			error(['Unexpected input to GetInputBoundConstraints(): ' varargin{argument_index} ])
		end
	end


end 