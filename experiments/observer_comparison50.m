function [results] = observer_comparison50(varargin)
	%observer_comparison50.m
	%
	%Description:
	%	The objective of this is to display the "movie" for leader follower dynamics.
	%
	%Assumptions:

	%%%%%%%%%%%%%%%%%%%
	%% Manage Inputs %%
	%%%%%%%%%%%%%%%%%%%

	switch nargin
	case 0
		1;
	case 2
		cube_x = varargin{1};
		cube_y = varargin{2};
	case 3
		cube_x = varargin{1};
		cube_y = varargin{2};
		fig_switch = varargin{3};
	case 4
		cube_x = varargin{1};
		cube_y = varargin{2};
		fig_switch = varargin{3};
		create_filters_flag = varargin{4};
	case 5
		cube_x = varargin{1};
		cube_y = varargin{2};
		fig_switch = varargin{3};
		create_filters_flag = varargin{4};
		debug_flag = varargin{5};
	otherwise
			error('Unexpected number of inputs.')
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	eta_w = 1.5; eta_v = 0.5;
	dt = 0.1;
	if ~exist('cube_x')
		cube_x = 2; cube_y = 2;
	end
	[~,ad0,foll_offsets] = get_lead_follow_aff_dyn3(cube_x,cube_y,dt, ...
													'disturb_info',eta_w,eta_v);

	n = size(ad0.A,1);

	Pu = 5*Polyhedron('lb',-ones(1,size(ad0.B,2)),'ub',ones(1,size(ad0.B,2)));

	T = 6;
	L1 = {ones(1,T)};
	%L2 = {[1,1,1],[0,1,1],[1,0,1],[1,1,0]}; periods_in_tunnel = 16; %Number of time periods in the tunnel.
	n0 = 5;
	L2 = {[1,1,1,ones(1,n0)],[0,1,1,ones(1,n0)],[0,0,1,ones(1,n0)],[1,0,1,ones(1,n0)],[0,1,0,ones(1,n0)],[1,1,0,ones(1,n0)]}; periods_in_tunnel = 5; %Number of time periods in the tunnel.
	T2 = length(L2{1});
	L3 = L1;
	T3 = 6;

	escape_periods = 2;
	settle_periods = 3;

	%Number of times to repeat frames
	rep_frames = 4;

	results.params.ad = ad0;
	results.params.Pu = Pu;
	results.params.L1 = L1;
	results.params.cube_x = cube_x;
	results.params.cube_y = cube_y;

	fs = 20; %Font Size
	ms = 118; %Marker Size. Default is 36
	follower_marker = 'h'; 	% Other optinos: 'd' or 'diamond'
							% 				 'p' or 'pentagram'
							%				 's' or 'square'
							%				 'h' or 'hexagram'

	nahs_results_image_folder = 'results/nahs2019/images/';

	if ~exist('debug_flag')
		debug_flag = 1;
	end

	%%%%%%%%%%%%%%%%%%%%
	%% Create Filters %%
	%%%%%%%%%%%%%%%%%%%%
	M2_temp = 4;
	M_t = M2_temp;
    
	if ~exist('create_filters_flag')
		try
			load(['results/nahs2019/lead_foll_startup_gains_' num2str(cube_x) 'by' num2str(cube_y) '.mat'])
			create_filters_flag = false;
		catch
			create_filters_flag = true;
		end
	end


	if create_filters_flag
		history.oo1 = {}; history.c1 = {};
		history.oo2 = {}; history.c2 = {};
		history.oo3 = {}; history.c3 = {};
		rec_num = 0;
		while M2_temp >= 1
			%Increment rec_num
			rec_num = rec_num + 1;

			[history.oo1{rec_num},history.c1{rec_num}] = ad0.rec_synthesis(	'Equalized','prefix','Minimize M2',M2_temp,L1, ...
																			'Pu',Pu, ...
																			'debug',debug_flag);
			[history.oo2{rec_num},history.c2{rec_num}] = ad0.rec_synthesis(	'Free','prefix','Minimize M3',M2_temp,history.oo1{rec_num}.M2,L1, ...
																			'Pu',Pu, ...
																			'debug',debug_flag);

			%Update M2
			M_t(1+(rec_num-1)*T+1:1+rec_num*T-1) = history.oo1{rec_num}.M2*ones(T-1,1);
			M_t(1+rec_num*T) = history.oo2{rec_num}.M3;

			%Update Loop Variables
			M2_temp = history.oo2{rec_num}.M3;
		end
		rec_tries = rec_num;

		M2_star = 1.3;
		[history.oo3,history.c3] = ad0.rec_synthesis('Equalized','prefix','Minimize M2',1,L2,'Pu',Pu);
		for period_idx = 1:periods_in_tunnel
			M_t(1+rec_tries*T+(period_idx-1)*T2+[1:T2-1]) = history.oo3.M2;
			M_t(1+rec_tries*T+period_idx*T2) = 1;
		end

		[history.oo_Lstar,history.c_Lstar] = ad0.rec_synthesis('Equalized','time','Minimize M2',1,L2,'Pu',Pu);

		[history.oo4,history.c4] = ad0.rec_synthesis('Equalized','prefix','Minimize M2',1,L3,'Pu',Pu);
		for period_idx = 1:escape_periods+settle_periods
			M_t(1+rec_tries*T+periods_in_tunnel*T2+(period_idx-1)*T3+[1:T3-1]) = history.oo4.M2;
			M_t(1+rec_tries*T+periods_in_tunnel*T2+period_idx*T3) = history.oo4.M2;
		end

		save(['results/nahs2019/lead_foll_startup_gains_' num2str(cube_x) 'by' num2str(cube_y) '.mat'],'history','M_t','M2_temp','rec_tries');
	else
		load(['results/nahs2019/lead_foll_startup_gains_' num2str(cube_x) 'by' num2str(cube_y) '.mat'])
	end

	%%%%%%%%%%%%%%%%%
	%% Create Path %%
	%%%%%%%%%%%%%%%%%

	l_0 = [-6;0];
	l_f = [25,0];
	l_targ = l_f + [10,6];

	l_x = [	linspace(l_0(1),0,T*rec_tries+1), linspace(0,l_f(1),periods_in_tunnel*T2), ...
				linspace(l_f(1),l_f(1)+1+2*cube_x+1,T3*escape_periods),linspace(l_f(1)+1+2*cube_x+1,l_targ(1),T3*settle_periods);
			linspace(l_0(2),0,T*rec_tries+1), linspace(0,l_f(2),periods_in_tunnel*T2), ...
				linspace(l_f(2),l_f(2),T3*escape_periods)			,linspace(l_f(2),l_targ(2),T3*settle_periods)];

	w = l_x(:,2:end) - l_x(:,1:end-1);

	

	results.history = history;
	results.params.M_Target = M2_temp;
	results.params.M_t = M_t;
	results.params.l_x = l_x;
	results.params.w_lead = w;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Plot all of the Agents %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	t_0 = 0;

	x0 = unifrnd(-M_t(1),M_t(1),n,1);
	x = x0;
	for contr_num = 1:rec_tries
		temp_traj = history.c2{contr_num}.simulate_1run( ad0 , M_t(1+contr_num*(T-1)) , 'in_sigma' , L1{1} , ...
																						'in_w' , w(:,(contr_num-1)*T+1:contr_num*T) , ...
																						'in_x0' , x(:,end)  );
		x = [x, temp_traj(:,2:end)];
	end

	for in_channel_idx = 1:periods_in_tunnel
		temp_traj = history.c3.simulate_1run( ad0 , M_t(1+contr_num*(T-1)) , 'in_w' , w(:,rec_tries*T+(in_channel_idx-1)*T2+1:rec_tries*T+(in_channel_idx)*T2) , ...
																			'in_x0' , x(:,end)  );
		x = [x, temp_traj(:,2:end)];
	end

	for period_idx = 1:escape_periods+settle_periods
		temp_traj = history.c4.simulate_1run(ad0, M_t(1+rec_tries*T+periods_in_tunnel*T2+(period_idx-1)*T3) , ...
											'in_w' , w(:,rec_tries*T+periods_in_tunnel*T2+(period_idx-1)*T3+[1:T3]) , ...
											'in_x0' , x(:,end)  );
		x = [x,temp_traj(:,2:end)];
	end

	results.x_start = x;

	foll_offsets2 = [foll_offsets(1,:)';foll_offsets(2,:)'];
	target_pos = repmat(foll_offsets2,1,size(l_x,2)) + [repmat(l_x(1,:),cube_x*cube_y,1); repmat(l_x(2,:),cube_x*cube_y,1) ];
	results.targets = target_pos;

	figure;
	hold on;
	for robot_idx = 1:cube_x*cube_y
		scatter(target_pos(robot_idx,:),target_pos(cube_x*cube_y+robot_idx,:))
	end

	x_true = x + target_pos;

	figure;
	for robot_idx = 1:cube_x*cube_y
		subplot(cube_x,cube_y,robot_idx)
		scatter(x_true(robot_idx,:),x_true(cube_x*cube_y+robot_idx,:))
	end

	t0 = 1;
	t_int = floor(size(x_true,2)/3);
	t_final = size(x_true,2);

	if ~exist('fig_switch')
		fig_switch = 5;
	end
	switch fig_switch
	case 1
		t_instances = [	t0;
						T*rec_tries+1;
						T*rec_tries+T2*10+1;
						T*rec_tries+T2*periods_in_tunnel+T3;
						t_final];
	case 2
		t_instances = [	t0;
						T*rec_tries+1;
						T*rec_tries+T2*10+1];
	case 3
		t_instances = [	T*rec_tries+T2*10+1;
						T*rec_tries+T2*periods_in_tunnel;
						T*rec_tries+T2*periods_in_tunnel+T3*escape_periods];
	case 4
		t_instances = [	T*rec_tries+T2*periods_in_tunnel+T3*escape_periods;
						T*rec_tries+T2*periods_in_tunnel+T3*(escape_periods+1);
						t_final];
	case 5
		t_instances = [ T*rec_tries+1;
						T*rec_tries+T2*floor(periods_in_tunnel/2);
						1+T*rec_tries+T2*periods_in_tunnel];
	case 6
		disp('Writing video... Need all time indices.')
	otherwise
		error('Unexpected figure switch input.')
	end

	%Plot Hyperboxes around the area where the followers can be.
	axis_lims = [l_0(1)-2*cube_x-6,l_targ(1)+2,-(l_targ(2)+2),l_targ(2)+2];

	if fig_switch <= 5
		figure('DefaultAxesFontSize',fs);
		for t_idx = 1:length(t_instances)
			subplot(1,length(t_instances),t_idx)
			hold on;
	        scatter(l_x(1,t_instances(t_idx)),l_x(2,t_instances(t_idx)),ms,'x')
			for robot_idx = 1:cube_x*cube_y
				scatter(x_true(robot_idx,t_instances(t_idx)), ...
						x_true(cube_x*cube_y+robot_idx,t_instances(t_idx)))
				axis([l_0(1)-6,l_f(1)+3,-10,10])
			end
		end

		figure('DefaultAxesFontSize',fs);
		for t_idx = 1:length(t_instances)
			subplot(1,length(t_instances),t_idx)
			hold on;
	        scatter(l_x(1,t_instances(t_idx)),l_x(2,t_instances(t_idx)),ms,'x')
			for robot_idx = 1:cube_x*cube_y
				scatter(x_true(robot_idx,t_instances(t_idx)), ...
						x_true(cube_x*cube_y+robot_idx,t_instances(t_idx)), ...
						ms,follower_marker)
				axis(axis_lims)
			end
			%Draw Error Box
			center_pos = l_x(:,t_instances(t_idx)) - [cube_x+1;0];
			box_width = 2*((cube_x-1) + M_t(t_instances(t_idx)));
			box_height = 2*(1+M_t(t_instances(t_idx)));
			rectangle('Position',[center_pos(1)-(1/2)*box_width,center_pos(2)-(1/2)*box_height, ...
									box_width, ...
									box_height ],...
						'FaceColor','none');
			%Create Obstacles
			rectangle('Position',[0,3,l_f(1)-0,axis_lims(4)-3], ...
						'FaceColor','red')
			rectangle('Position',[0,axis_lims(3),l_f(1)-0,-3-axis_lims(3)], ...
						'FaceColor','red')

		end
		set(gcf,'Units','Normalized','Position',[0,0,1,1])
		saveas(gcf,[ nahs_results_image_folder 'consensus_following_' num2str(cube_x) 'by' num2str(cube_y) '_v' num2str(fig_switch) ], 'epsc')

		if fig_switch == 5
			figure('DefaultAxesFontSize',fs);
			for t_idx = 1:length(t_instances)-1
				subplot(1,length(t_instances)-1,t_idx)
				hold on;
		        scatter(l_x(1,t_instances(t_idx)),l_x(2,t_instances(t_idx)),ms,'x')
				for robot_idx = 1:cube_x*cube_y
					scatter(x_true(robot_idx,t_instances(t_idx)), ...
							x_true(cube_x*cube_y+robot_idx,t_instances(t_idx)), ...
							ms, follower_marker)
					axis(axis_lims)
				end
				%Draw Error Box
				center_pos = l_x(:,t_instances(t_idx)) - [cube_x+1;0];
				if t_idx == 1
					box_width = 2*((cube_x-1) + M_t(t_instances(t_idx)));
					box_height = 2*(1+M_t(t_instances(t_idx)));
				else
					box_width = 2*((cube_x-1) + M_t(t_instances(t_idx)));
					box_height = 2*(1+history.oo_Lstar.M2 );
				end

				%Create Obstacles
				rectangle('Position',[0,2.5,l_f(1)-0,axis_lims(4)-3], ...
							'FaceColor','red')
				rectangle('Position',[0,axis_lims(3),l_f(1)-0,-2.5-axis_lims(3)], ...
							'FaceColor','red')
				%Create bounding box
				rectangle('Position',[center_pos(1)-(1/2)*box_width,center_pos(2)-(1/2)*box_height, ...
									box_width, ...
									box_height ],...
						'FaceColor','none');
				set(gcf,'Units','Normalized','Position',[0,0,1,1])
				saveas(gcf,[ nahs_results_image_folder 'consensus_following_' num2str(cube_x) 'by' num2str(cube_y) '_worst_case'], 'epsc')
			end
		end
	end

	%% Movie Writing
	if fig_switch == 6
		%Create figure
		figure;

		%Create video writer
		v = VideoWriter([ nahs_results_image_folder 'formation_vid0.mp4'],'MPEG-4');
		open(v);

		for t_idx = 1:size(x_true,2)
			%Plot leader
			scatter(l_x(1,t_idx),l_x(2,t_idx),ms,'x')
			hold on;
			
			%Plot followers
			for robot_idx = 1:cube_x*cube_y
				scatter(x_true(robot_idx,t_idx), ...
						x_true(cube_x*cube_y+robot_idx,t_idx), ...
						ms,follower_marker)
				axis(axis_lims)
			end
			%Draw Error Box
			center_pos = l_x(:,t_idx) - [cube_x+1;0];
			
			box_width = 2*((cube_x-1) + M_t(t_idx));
			box_height = 2*(1+M_t(t_idx));
			

			%Create Obstacles
			rectangle('Position',[0,2.5,l_f(1)-0,axis_lims(4)-3], ...
						'FaceColor','red')
			rectangle('Position',[0,axis_lims(3),l_f(1)-0,-2.5-axis_lims(3)], ...
						'FaceColor','red')
			%Create bounding box
			rectangle('Position',[center_pos(1)-(1/2)*box_width,center_pos(2)-(1/2)*box_height, ...
								box_width, ...
								box_height ],...
					'FaceColor','none');
			set(gcf,'Units','Normalized','Position',[0,0,1,1])

			%Write Video Frame
			frame = getframe(gcf);
			for rep_idx = 1:rep_frames
   				writeVideo(v,frame);
   			end
   			hold off;
		end
		
		close(v);

	end


end