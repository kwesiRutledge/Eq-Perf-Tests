function [ Q , r , Pi1 , Piu ] = get_rec_optim_vars(varargin)
	%Description:
	%	Gets the variables necessary to define a finite horizon reachability problem.
	%	- Q-Parameterized versions of the feedback gains (F,f) [Notation used is the NAHS submission variety.]
	%	- Dual variable that satisfy polytope inclusion of the state
	%	- Dual variable that is used to define satisfaction of the input constraint.
	%
	%Usage:
	%	[ Q , r , Pi1 , Piu ] = cg.get_rec_optim_vars(lcsas,word,P_x0)
	%
	%Inputs:
	%	lcsas - A language constrained switched affine system (LCSAS) defined in an LCSAS object.
	%	word - A word that is feasible/allowed given the definition of the lcsas.

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	cg = varargin{1};

	if isa(varargin{2},'LCSAS')
		sys = varargin{2};
	elseif isa(varargin{2},'Aff_Dyn')
		error('This function currently does not support Aff_Dyn objects. Maybe in the future?')
	else
		error('Unexpected input for the second input. Expected an LCSAS or Aff_Dyn type of object.')
	end

	word = varargin{3};
	P_x0 = varargin{4};

	if ~isa(P_x0,'Polyhedron')
		error('P_x0 must be given as a Polyhedron.')
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	n = size(sys.Dyn(1).A,1);
	m = size(sys.Dyn(1).B,2);
	p = size(sys.Dyn(1).C,1);
	wd = size(sys.Dyn(1).B_w,2);
	vd = size(sys.Dyn(1).C_v,2);

	T_i = length(word);

	%%%%%%%%%%%%%%%%%
	%% Define Q,u0 %%
	%%%%%%%%%%%%%%%%%

	Q = sdpvar(m*T_i,p*T_i,'full');
	r = sdpvar(m*T_i,1,'full');

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Define the Dual Variables %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	Pi1 = sdpvar(size(P_x0.A,1),2*(wd+vd)*T_i+2*n,'full');

end