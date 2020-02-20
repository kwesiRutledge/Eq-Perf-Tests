
function [results] = observer_comparison56( varargin )
	%Description:
	%	Exemplifying the 'constraints' on distinguishability.

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	switch nargin
		case 1
			plot_traj2 = varargin{1};
	end

	if ~exist('plot_traj2')
		plot_traj2 = false;
	end

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

	eta_v = 0; eta_w = 0.0;
	Pv1 = Polyhedron('lb',-eta_v*ones(1,dim) ,'ub',eta_v*ones(1,dim));
	Pw1 = Polyhedron('lb',-eta_w*ones(1,dim) ,'ub',eta_w*ones(1,dim));

	ad1 = Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1);

	f2 = [1;0];
	f3 = -f1;
	f4 = -f2;

	lcss1 = [	ad1, ...
				Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1), ...
				Aff_Dyn(A1,B1,f3,C1,Pw1,Pv1), ...
				Aff_Dyn(A1,B1,f4,C1,Pw1,Pv1)];

	results.exp1.L = L1;
	results.exp1.lcsas = lcss1;

	%Define Sets
	eta_u = 0; eta_x0 = 0.3;
	P_u = Polyhedron('lb',-eta_u*ones(1,dim) ,'ub',eta_u*ones(1,dim));
	P_x0 = Polyhedron('lb',-eta_x0*ones(1,dim),'ub',eta_x0*ones(1,dim));

	% Select Matrix
	RT = @(n,t,T) [zeros(n,n*t) eye(n) zeros(n,n*(T-t))];

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Plotting the futures of the 3 Languages %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	
	timer_start = tic;

	Hc = {}; Sc = {}; Jc = {}; fc = {};

	for word_idx = 1:length(L1.words)
		[Hc{word_idx},Sc{word_idx},~,Jc{word_idx},fc{word_idx}] = get_mpc_matrices(lcss1,'word',L1.words{word_idx});
		%Create Cartesian Product of Polyhedra
		P_wT{word_idx} = 1; P_vT{word_idx} = 1; P_uT{word_idx} = 1;
		for sym_idx = 1:length(L1.words{word_idx})
			P_wT{word_idx} = P_wT{word_idx} * lcss1(L1.words{word_idx}(sym_idx)).P_w;
			P_vT{word_idx} = P_vT{word_idx} * lcss1(L1.words{word_idx}(sym_idx)).P_v;
			P_uT{word_idx} = P_uT{word_idx} * P_u;
		end
	end

	T = L1.find_longest_length();

	%% Create P_wT and P_vT
	figure;
	for seq_ind = 1:length(Hc)
		X = (Hc{seq_ind}*(P_wT{seq_ind}+fc{seq_ind})) + Sc{seq_ind}*P_uT{seq_ind} + Jc{seq_ind}*P_x0;
		subplot(1,length(Hc),seq_ind)
		hold on;
		plot(RT(dim,0,L1.find_longest_length())*X,'Color','white')
		plot(RT(dim,1,L1.find_longest_length())*X,'Color','cyan')
		plot(RT(dim,2,L1.find_longest_length())*X,'Color','magenta')
		plot(RT(dim,3,L1.find_longest_length())*X,'Color','orange')
		xlabel('$x_1(t)$','Interpreter','latex')
		ylabel('$x_2(t)$','Interpreter','latex')
		axis([-2 3 -2 3])
	end

	timing_info.plot1 = toc(timer_start);

	saveas(gcf,'results/oc56_trajectories_sys1.png')

	%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Create Belief Graph %%
	%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('- Creating BeliefGraph 1.')

	timer_start = tic;
	bg = BeliefGraph(lcss1,L1, P_u, P_x0, 'fb_method','state');
	timing_info.create_bg1 = toc(timer_start);

	figure;
	bg.plot();
	disp('- Plotted BeliefGraph 1.')

	saveas(gcf,'results/oc56_belief_graph_sys1.png')

	results.exp1.bg = bg;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Create Belief Graph with More Edges %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	clear bg eta_v eta_w ad1 Pv1 Pw1 P_u P_wT P_vT P_x0;

	eta_v = 0.2; eta_w = 0.2;
	Pv1 = Polyhedron('lb',-eta_v*ones(1,dim) ,'ub',eta_v*ones(1,dim));
	Pw1 = Polyhedron('lb',-eta_w*ones(1,dim) ,'ub',eta_w*ones(1,dim));

	eta_u = 10^(-3); eta_x0 = 0.1;
	P_u = Polyhedron('lb',-eta_u*ones(1,dim) ,'ub',eta_u*ones(1,dim));
	P_x0 = Polyhedron('lb',-eta_x0*ones(1,dim),'ub',eta_x0*ones(1,dim));


	ad1 = Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1);

	f2 = [1;0];
	f3 = -f1;
	f4 = -f2;

	lcss2 = [	ad1,...
				Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1),...
				Aff_Dyn(A1,B1,f3,C1,Pw1,Pv1),...
				Aff_Dyn(A1,B1,f4,C1,Pw1,Pv1)];

	disp('- Created Switching Aff_Dyn system.')

	timer_start = tic; %Timing the plotting?

	Hc = {}; Sc = {}; Jc = {}; fc = {};

	for word_idx = 1:length(L1.words)
		[Hc{word_idx},Sc{word_idx},Cc{word_idx},Jc{word_idx},fc{word_idx}] = get_mpc_matrices(lcss2,'word',L1.words{word_idx});

		%Create Cartesian Product of Polyhedra
		P_wT{word_idx} = 1; P_vT{word_idx} = 1; P_uT{word_idx} = 1;
		for sym_idx = 1:length(L1.words{word_idx})
			P_wT{word_idx} = P_wT{word_idx} * lcss2(L1.words{word_idx}(sym_idx)).P_w;
			P_vT{word_idx} = P_vT{word_idx} * lcss2(L1.words{word_idx}(sym_idx)).P_v;
			P_uT{word_idx} = P_uT{word_idx} * P_u;
		end
	end

	T = L1.find_longest_length();

	if plot_traj2
		figure;
		for seq_ind = 1:L1.cardinality()
			
			X = (Hc{seq_ind}*(P_wT{seq_ind}+fc{seq_ind})) + Jc{seq_ind}*P_x0;
			X.computeVRep;
			Y = Cc{seq_ind}*X + P_vT{seq_ind};

			subplot(1,length(Hc),seq_ind)
			hold on;
			plot(RT(dim,0,T-1)*Y,'Color','white')
			plot(RT(dim,1,T-1)*Y,'Color','cyan')
			plot(RT(dim,2,T-1)*Y,'Color','magenta')
			plot(RT(dim,3,T-1)*Y,'Color','orange')
			%title(['Trajectory of Word ' num2str(seq_ind)])
			xlabel('$x_1(t)$','Interpreter','latex')
			ylabel('$x_2(t)$','Interpreter','latex')
			axis([-2 3 -2 3])
		end
	end

	timing_info.plot2 = toc(timer_start); %Timing	

	% saveas(gcf,'results/oc56_trajectories_sys2.png')
	disp('- Creating BeliefGraph 2.')

	timer_start = tic;
	bg = BeliefGraph(lcss2,L1, P_u, P_x0,'fb_method','output','accel_flag',true);
	timing_info.create_bg2 = toc(timer_start);

	figure;
	bg.plot();
	disp('Plotted BeliefGraph 2.')

	saveas(gcf,'results/oc56_belief_graph_sys2.png')

	results.exp2.bg2 = bg;

	results.timing_info = timing_info;

end
