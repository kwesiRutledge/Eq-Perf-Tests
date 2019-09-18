function [results] = observer_comparison57( varargin )
	%Description:
	%	Exemplifying the 'constraints' on distinguishability.

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	dim = 2;

	L1 = Language([1,1,2,3],[1,2,2,4],[1,2,4,3]);

	%Create a simple Language Constrainted Switching System
	A1 = eye(dim);
	B1 = eye(dim);
	C1 = eye(dim);
	f1 = [0;1];

	eta_v = 0; eta_w = 0.5;
	Pv1 = Polyhedron('lb',-eta_v*ones(1,dim) ,'ub',eta_v*ones(1,dim));
	Pw1 = Polyhedron('lb',-eta_w*ones(1,dim) ,'ub',eta_w*ones(1,dim));

	eta_u = 0.1;
	P_u = Polyhedron('lb',-eta_u*ones(1,dim) ,'ub',eta_u*ones(1,dim));

	ad1 = Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1);

	f2 = [1;0];
	f3 = -f1;
	f4 = -f2;

	lcss = [	ad1,...
				Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1),...
				Aff_Dyn(A1,B1,f3,C1,Pw1,Pv1),...
				Aff_Dyn(A1,B1,f4,C1,Pw1,Pv1)];

	%Define Sets
	eta_u = 0; eta_x0 = 0.3;
	P_u = Polyhedron('lb',-eta_u*ones(1,dim) ,'ub',eta_u*ones(1,dim));
	P_x0 = Polyhedron('lb',-eta_x0*ones(1,dim),'ub',eta_x0*ones(1,dim));

	% Select Matrix
	RT = @(n,t,T) [zeros(n,n*t) eye(n) zeros(n,n*(T-t))];

	%% Synthesizing a feedback gain that tries to recover %%

	disp('Experiment 1: Generating constraints')
	lcss3 = LCSAS(lcss,L1);
	disp('- Created LCSAS object.')
	% bg3 = BeliefGraph(lcss2,L1, P_u, P_x0);
	% blang3 = bg3.get_belief_language();
	[contr3, opt_out3, constrs] = lcss3.rec_synthesis(P_x0, 'P_u', P_u);
	disp('- Finished Recovery Synthesis.')

	results.system = lcss3;
	results.controller = contr3;
	results.optim_out = opt_out3;
	results.constraints = constrs;

end