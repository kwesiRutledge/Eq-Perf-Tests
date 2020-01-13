function [results] = observer_comparison60( varargin )
	%observer_comparison60.m
	%Description:
	%	Testing various systems with the algorithms that have been developed so far.
	%

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if nargin >= 1
		load_data_flag = varargin{1};
	end

	if nargin >= 3
		c_sq.dim_x = varargin{2};
		c_sq.dim_y = varargin{3};
	end


	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	if ~exist('c_sq')
		c_sq.dim_x = 2;
		c_sq.dim_y = 2;
	end

	dt = 0.1;

	L = Language([1,2,2,2],[1,1,2,2],[1,1,1,2]);
	in_sys = get_consensus_dyn(c_sq.dim_x,c_sq.dim_y,dt,'L',L );

	n_x = size(in_sys.Dyn(1).A,1);
	n_u = size(in_sys.Dyn(1).B,2);

	eta_u = 100.0;
	Pu = Polyhedron('lb',-eta_u*ones(1,n_u),'ub',eta_u*ones(1,n_u));

	eta_x = 2.0;
	Px0 = Polyhedron('lb',-eta_x*ones(1,n_x),'ub',eta_x*ones(1,n_x));	

	results.params.sys = in_sys;
	results.params.Pu = Pu;
	results.params.Px0 = Px0;

	verbosity = 0; %Verbosity of Functions. Gives debugging info

	%Data file parameters.
	save_file_name = ['data/oc60_interm_results_' num2str(c_sq.dim_x) 'x' num2str(c_sq.dim_y) '.mat'];
	if ~exist('load_data_flag')
		load_data_flag = false;
	end
	run_opt_flag = true;

	if load_data_flag
		load(save_file_name)
		load_data_flag = true;
	end

	ops = sdpsettings('verbose',verbosity);

	results_filename_base = 'results/oc60_results_';

	%%%%%%%%%%%
	%% Tests %%
	%%%%%%%%%%%

	if (~load_data_flag) && run_opt_flag
		[BG,contr,opt_data] = in_sys.synth_robust_reach_contr(Px0,Pu,2*Px0, 'debug' , verbosity)
		save(save_file_name,'contr','opt_data','BG','c_sq','eta_u','eta_x')
	elseif load_data_flag && run_opt_flag
		[~,contr,opt_data] = in_sys.synth_robust_reach_contr(Px0,Pu,1.5*Px0,'BG',BG, 'debug' , verbosity)
		save(save_file_name,'contr','opt_data','BG','c_sq','eta_u','eta_x')
	end

	results.BG = BG;
	results.contr = contr;
	results.opt_data = opt_data;

	figure;
	BG.plot();
	saveas(gcf,[results_filename_base 'belief_graph'], 'epsc')


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Simulate the System With this Form of Controller %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	in_w = [ unifrnd(-0.5,0.5,2,2), 0.75*ones(2,2) ];
	sig = [1,1,2,2];

	in_x0 = [-2;2;2;-2];

	[x, u , y , ~ ] = contr.simulate_1run(in_sys,Px0, 'in_sigma' , sig , 'in_w' , in_w , 'in_x0' , in_x0 );
	max_x = max(max(x([1:2:end],:)));
	min_x = max(max(-x([1:2:end],:)));
	min_x = -min_x;

	max_y = max(max(x([2:2:end],:)));
	min_y = max(max(-x([2:2:end],:)));
	min_y = -min_y;

	%Rectangle Plotting Info
	center_pos = [0;0];
	M2 = opt_data.M2;
	M3 = opt_data.M3;

	M2_prime = M2 + 0.3;

	figure;
	hold on;
	plot_idx = 0;
	scatter(x(1,:),x(2,:))
	rectangle('Position',[center_pos(1)-M2,center_pos(2)-M2, ...
									2*M2, ...
									2*M2 ],...
						'FaceColor','none');
	rectangle('Position',[center_pos(1)-M3,center_pos(2)-M3, ...
									2*M3, ...
									2*M3 ],...
						'FaceColor',[0.3010 0.7450 0.9330],'LineStyle', '--')

	axis([-M2_prime, M2_prime, -M2_prime, M2_prime])
	%title(['$t=' num2str(t-1) '$'],'Interpreter','latex')

	saveas(gcf,[results_filename_base 'comic_strip'], 'epsc')

	results.x = x;
	results.u = u;
	results.y = y;
	results.sig = sig;
	results.b = contr.b_hist

	% figure;
	% plot3(x(1,:),x(2,:),x(3,:))

	% figure;
	% plot(x(:,1),x(:,2))

	% T = size(x,2);
	% figure;
	% hold on;
	% plot([1:T],x(1,:))
	% plot([1:T],x(2,:))
	% plot([1:T],x(3,:))

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Save The Larger Variables to a Data File %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end