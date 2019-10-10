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

	verbosity = 1; %Verbosity of Functions. Gives debugging info

	%Data file parameters.
	save_file_name = ['data/oc60_interm_results_ ' num2str(c_sq.dim_x) 'x' num2str(c_sq.dim_y) '.mat'];
	if ~exist('load_data_flag')
		load_data_flag = false;
	end
	run_opt_flag = false;

	if load_data_flag
		load(save_file_name)
		load_data_flag = true;
	end

	cg = constr_gen();

	ops = sdpsettings('verbose',verbosity);

	%%%%%%%%%%%
	%% Tests %%
	%%%%%%%%%%%

	if ~load_data_flag
		[BG,contr,opt_data] = in_sys.synth_robust_reach_contr(Px0,Pu,2*Px0)
		save(save_file_name,'contr','opt_data','BG','c_sq')
	end

	results.BG = BG;
	results.contr = contr;
	results.opt_data = opt_data;


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Simulate the System With this Form of Controller %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	in_w = [ unifrnd(-0.5,0.5,2,4), ones(2,2) ];
	sig = [1,1,1,1,2,2];

	[x, u , y , sig] = contr.simulate_1run(in_sys,Px0, 'in_sigma' , sig , 'in_w' , in_w );
	figure;
	plot_idx = 0;
	for t = [0,3,6]
		plot_idx = plot_idx + 1;
		subplot(1,3,plot_idx)
		hold on;
		for agent_idx = 1:c_sq.dim_x*c_sq.dim_y
			scatter(x((agent_idx-1)*2+1,t+1),x((agent_idx-1)*2+2,t+1))
		end
	end
	

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