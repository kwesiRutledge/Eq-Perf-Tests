function [results] = observer_comparison56( varargin )
	%Description:
	%	Exemplifying the 'constraints' on distinguishability.

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	dim = 2;

	L1 = Language([1,1,2],[1,2,2],[1,2,4]);

	%Create a simple Language Constrainted Switching System
	A1 = eye(dim);
	B1 = eye(dim);
	C1 = eye(dim);
	f1 = [0;1];

	eta_v = 0; eta_w = 0.02;
	Pv1 = Polyhedron('lb',-eta_v*ones(1,dim) ,'ub',eta_v*ones(1,dim));
	Pw1 = Polyhedron('lb',-eta_w*ones(1,dim) ,'ub',eta_w*ones(1,dim));

	ad1 = Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1);

	f2 = [1;0];
	f3 = -f1;
	f4 = -f2;

	lcss1 = [	ad1,...
				Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1),...
				Aff_Dyn(A1,B1,f3,C1,Pw1,Pv1),...
				Aff_Dyn(A1,B1,f4,C1,Pw1,Pv1)];

	results.parameters.L1 = L1;
	results.parameters.lcsas = lcss1;

	%Define Sets
	eta_u = 0; eta_x0 = 0.3;
	P_u = Polyhedron('lb',-eta_u*ones(1,dim) ,'ub',eta_u*ones(1,dim));
	P_x0 = Polyhedron('lb',-eta_x0*ones(1,dim),'ub',eta_x0*ones(1,dim));

	% Select Matrix
	RT = @(n,t,T) [zeros(n,n*t) eye(n) zeros(n,n*(T-t))];

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Plotting the futures of the 3 Languages %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	Hc = {}; Sc = {}; Jc = {}; fc = {};

	for word_idx = 1:length(L1.words)
		[Hc{word_idx},Sc{word_idx},~,Jc{word_idx},fc{word_idx}] = get_mpc_matrices(lcss1,'word',L1.words{word_idx});
	end

	T = L1.find_longest_length();

	figure;
	for seq_ind = 1:length(Hc)
		X = (Hc{seq_ind}*(Pw1*Pw1*Pw1+fc{seq_ind})) + Sc{seq_ind}*(P_u*P_u*P_u) + Jc{seq_ind}*P_x0;
		subplot(1,length(Hc),seq_ind)
		hold on;
		plot(RT(dim,0,3)*X,'Color','white')
		plot(RT(dim,1,3)*X,'Color','cyan')
		plot(RT(dim,2,3)*X,'Color','magenta')
		plot(RT(dim,3,3)*X,'Color','orange')
		axis([-5 5 -5 5])
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Create Belief Graph %%
	%%%%%%%%%%%%%%%%%%%%%%%%%

	bg = BeliefGraph(lcss1,L1, P_u, P_x0);

	figure;
	bg.plot();

	results.bg1 = bg;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Create Belief Graph with More Edges %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	clear bg eta_v eta_w

	eta_v = 0; eta_w = 0.4;
	Pv1 = Polyhedron('lb',-eta_v*ones(1,dim) ,'ub',eta_v*ones(1,dim));
	Pw1 = Polyhedron('lb',-eta_w*ones(1,dim) ,'ub',eta_w*ones(1,dim));

	eta_u = 0.2;
	P_u = Polyhedron('lb',-eta_u*ones(1,dim) ,'ub',eta_u*ones(1,dim));

	ad1 = Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1);

	f2 = [1;0];
	f3 = -f1;
	f4 = -f2;

	lcss2 = [	ad1,...
				Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1),...
				Aff_Dyn(A1,B1,f3,C1,Pw1,Pv1),...
				Aff_Dyn(A1,B1,f4,C1,Pw1,Pv1)];

	Hc = {}; Sc = {}; Jc = {}; fc = {};

	for word_idx = 1:length(L1.words)
		[Hc{word_idx},Sc{word_idx},~,Jc{word_idx},fc{word_idx}] = get_mpc_matrices(lcss2,'word',L1.words{word_idx});
	end

	T = L1.find_longest_length();

	figure;
	for seq_ind = 1:length(Hc)
		X = (Hc{seq_ind}*(Pw1*Pw1*Pw1+fc{seq_ind})) + Sc{seq_ind}*(P_u*P_u*P_u) + Jc{seq_ind}*P_x0;
		subplot(1,length(Hc),seq_ind)
		hold on;
		plot(RT(dim,0,3)*X,'Color','white')
		plot(RT(dim,1,3)*X,'Color','cyan')
		plot(RT(dim,2,3)*X,'Color','magenta')
		plot(RT(dim,3,3)*X,'Color','orange')
		axis([-5 5 -5 5])
	end

	bg = BeliefGraph(lcss2,L1, P_u, P_x0);

	figure;
	bg.plot();

	results.bg2 = bg;

	%% Synthesizing a feedback gain that tries to place the user in a 

end
