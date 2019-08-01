function [results] = observer_comparison50(varargin)
	%observer_comparison50.m
	%
	%Description:
	%	The objective of this is to display the "movie" for leader follower dynamics.
	%
	%Assumptions:

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	eta_w = 0.4; eta_v = 0.1;
	dt = 0.1;
	cube_x = 2; cube_y = 2;
	[~,ad0,foll_offsets] = get_lead_follow_aff_dyn3(cube_x,cube_y,dt, ...
													'disturb_info',eta_w,eta_v);

	n = size(ad0.A,1);

	Pu = 5*Polyhedron('lb',-ones(1,size(ad0.B,2)),'ub',ones(1,size(ad0.B,2)));

	T = 6;
	L1 = {ones(1,T)};
	L2 = {[1,1,1],[0,1,1],[1,0,1],[1,1,0]};

	periods_in_tunnel = 3; %Number of time periods in the tunnel.

	rec_tries = 3;

	results.params.ad = ad0;
	results.params.Pu = Pu;
	results.params.L1 = L1;

	%%%%%%%%%%%%%%%%%
	%% Create Path %%
	%%%%%%%%%%%%%%%%%

	3 - 2;

	l0 = [6;0];

	l_x0 = repmat(l0,1,1+rec_tries*T+periods_in_tunnel*3)+ ...
						[ fliplr(linspace(-l0(1),4,1+rec_tries*T+periods_in_tunnel*3)) ; zeros(1,1+rec_tries*T+periods_in_tunnel*3)];

	w = l_x0(:,2:end) - l_x0(:,1:end-1);

	M2_temp = 4.8;
	M_t = M2_temp;

	create_filters_flag = false;

	if create_filters_flag
		history.oo1 = {}; history.c1 = {};
		history.oo2 = {}; history.c2 = {};
		history.oo3 = {}; history.c2 = {};
		for rec_num = 1:rec_tries
			[history.oo1{rec_num},history.c1{rec_num}] = ad0.rec_synthesis('Equalized','prefix','Minimize M2',M2_temp,L1,'Pu',Pu);
			[history.oo2{rec_num},history.c2{rec_num}] = ad0.rec_synthesis('Free','prefix','Minimize M3',M2_temp,history.oo1{rec_num}.M2,L1,'Pu',Pu);

			%Update M2
			M_t(1+(rec_num-1)*(T-1)+1:1+rec_num*(T-1)-1) = history.oo1{rec_num}.M2*ones(T-2,1);
			M_t(1+rec_num*(T-1)) = history.oo2{rec_num}.M3;

			M2_temp = history.oo2{rec_num}.M3;
		end

		[history.oo3,history.c3] = ad0.rec_synthesis('Equalized','prefix','Feasible Set',1,history.oo1{rec_num}.M2,L2,'Pu',Pu);

		save('results/nahs2019/lead_foll_startup_gains.mat','history','M_t','M2_temp');
	else
		load('results/nahs2019/lead_foll_startup_gains.mat')
	end

	results.params.M_Target = M2_temp;
	results.params.M_t = M_t;
	results.params.l_x = l_x0;
	results.params.w_lead = w;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Plot all of the Agents %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	t_0 = 0;

	x0 = unifrnd(-M_t(1),M_t(1),n,1);
	x = x0;
	for contr_num = 1:length(history.c2)
		x = [x, history.c2{contr_num}.simulate_1run( ad0 , M_t(1+contr_num*(T-1)) , 'in_sigma' , L1{1} , ...
																					'in_w' , w(:,(contr_num-1)*T+1:contr_num*T) , ...
																					'in_x0' , x(:,end)  )]
	end

	for in_channel_idx = 1:periods_in_tunnel
		x = [x, history.c3.simulate_1run( ad0 , M_t(1+contr_num*(T-1)) , 'in_w' , w(:,length(history.c2)*T+(in_channel_idx-1)*3+1:length(history.c2)*T+(in_channel_idx)*3) , ...
																		'in_x0' , x(:,end)  )]
	end

	results.x_start = x;

	foll_offsets2 = [foll_offsets(1,:)';foll_offsets(2,:)'];
	target_pos = repmat(foll_offsets2,1,size(l_x0,2)) + [repmat(l_x0(1,:),cube_x*cube_y,1); repmat(l_x0(2,:),cube_x*cube_y,1) ];
	results.targets = target_pos;

	figure;
	hold on;
	for robot_idx = 1:cube_x*cube_y
		scatter(target_pos(robot_idx,:),target_pos(cube_x*cube_y+robot_idx,:))
	end

end