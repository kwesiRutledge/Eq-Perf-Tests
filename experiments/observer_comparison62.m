function [results] = observer_comparison62( varargin )
	%observer_comparison62.m
	%Description:
	%	Testing the spacecraft lcsas control with our current algorithms.

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if nargin >= 1
		load_data_flag = varargin{1};
	end

	if nargin >= 2
		L_len = varargin{2};
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	if ~exist('L_len')
		L_len = 6;
	end

	%Data file parameters.
	save_file_name = ['data/oc62_interm_results_spacecraft_Llen' num2str(L_len) '.mat'];
	if ~exist('load_data_flag')
		load_data_flag = false;
	end
	run_opt_flag = true;

	if ~exist('c_sq')
		c_sq.dim_x = 2;
		c_sq.dim_y = 2;
	end

	if ~load_data_flag

		dt = 0.1;

		switch L_len
			case 6
				in_sys = get_spacecraft_lcsas1();
			case 5
				L = Language([1,2,1,1,1],[1,3,1,1,1],[1,2,3,1,1]);
				in_sys = get_spacecraft_lcsas1('L',L);
			otherwise
				error(['Unexpected length given: ' num2str(L_len)])
			
		end

		n_x = size(in_sys.Dyn(1).A,1);
		n_u = size(in_sys.Dyn(1).B,2);

		eta_u = 3.0;
		Pu = Polyhedron('lb',-eta_u*ones(1,n_u),'ub',eta_u*ones(1,n_u));

		eta_x = 0.5;
		Px0 = Polyhedron('lb',-eta_x*ones(1,n_x),'ub',eta_x*ones(1,n_x));		
	
	else
		load(save_file_name);

		in_sys = BG.lcsas;

		n_x = size(in_sys.Dyn(1).A,1);
		n_u = size(in_sys.Dyn(1).B,2);

		Pu = Polyhedron('lb',-eta_u*ones(1,n_u),'ub',eta_u*ones(1,n_u));
		Px0 = Polyhedron('lb',-eta_x*ones(1,n_x),'ub',eta_x*ones(1,n_x));	

	end

	verbosity = 1; %Verbosity of Functions. Gives debugging info
	ops = sdpsettings('verbose',verbosity);

	results.params.sys = in_sys;
	results.params.Pu = Pu;
	results.params.Px0 = Px0;

	%%%%%%%%%%%
	%% Tests %%
	%%%%%%%%%%%

	x0_offset = [1.5;1.5;1;-1.5;1.5;1];

	if (~load_data_flag) && run_opt_flag
		[BG,contr,opt_data] = in_sys.synth_robust_reach_contr(Px0,Pu)
		save(save_file_name,'contr','opt_data','BG','c_sq','eta_u','eta_x')
	elseif load_data_flag && run_opt_flag
		[~,contr,opt_data] = in_sys.synth_robust_reach_contr(Px0,Pu,'BG',BG)%,'P_target',(5/4)*Px0+x0_offset)
		save(save_file_name,'contr','opt_data','BG','c_sq','eta_u','eta_x')
	end

	results.BG = BG;
	results.contr = contr;
	results.opt_data = opt_data;


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Simulate the System With this Form of Controller %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%in_w = [ unifrnd(-0.2,0.2,2,5) ];
	in_w = [ unifrnd(-0.2,0.2,2,1), [0;0.3] ,unifrnd(-0.2,0.2,2,3) ];
	sig = [1,2,1,1,1];

	[x, u , y , sig] = contr.simulate_1run(in_sys,Px0, 'in_sigma' , sig , 'in_w' , in_w );
	figure; hold on;
	plot_idx = 0;
	scatter(x(1,:),x(4,:))
	

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