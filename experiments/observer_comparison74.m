function [results] = observer_comparison74( varargin )
	%observer_comparison74.m
	%Description:
	%	Plots the results from our drone tests using an image of a drone.
	%	For ECC 2020 Presentation

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	argin_idx = 1;
	while argin_idx <= nargin
		switch varargin{argin_idx}
			case 'load_data_flag'
				load_data_flag = varargin{argin_idx+1};
				argin_idx = argin_idx + 2;
			case 'c_sq'
				c_sq.dim_x = varargin{argin_idx+1};
				c_sq.dim_y = varargin{argin_idx+2};
				argin_idx = argin_idx + 3;
			case 'L'
				L = varargin{argin_idx+1};
				argin_idx = argin_idx + 2;
			case 'disturb_bounds'
				eta_w = varargin{argin_idx+1};
				eta_v = varargin{argin_idx+2};
				argin_idx = argin_idx + 3;
			case 'eta_u'
				eta_u = varargin{argin_idx+1};
				argin_idx = argin_idx + 2;
			case 'save_file_name'
				save_file_name = varargin{argin_idx + 1};
				argin_idx = argin_idx + 2;
			otherwise
				error(['Unexpected input to the experiment: ' varargin{argin_idx} ])
		end
	end


	if ~exist('c_sq')
		c_sq.dim_x = 2;
		c_sq.dim_y = 2;
	end

	dt = 0.1;

	if ~exist('L')
		L = Language([1,2,2,1,1],[1,1,2,2,1],[1,1,1,2,2]);
	end

	if ~exist('eta_v')
		eta_v = 0.2;
	end

	if ~exist('eta_w')
		eta_w = 0.2;
	end

	if ~exist('load_data_flag')
		load_data_flag = false;
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	in_sys = get_consensus_dyn(c_sq.dim_x,c_sq.dim_y,dt,'L',L,'disturb_params', eta_w , eta_v );

	n_x = size(in_sys.Dyn(1).A,1);
	n_u = size(in_sys.Dyn(1).B,2);

	if ~exist('eta_u')
		eta_u = 100.0;
	end
	Pu = Polyhedron('lb',-eta_u*ones(1,n_u),'ub',eta_u*ones(1,n_u));

	eta_x = 2.0;
	Px0 = Polyhedron('lb',-eta_x*ones(1,n_x),'ub',eta_x*ones(1,n_x));	

	results.params.sys = in_sys;
	results.params.Pu = Pu;
	results.params.Px0 = Px0;

	verbosity = 0; %Verbosity of Functions. Gives debugging info

	%Data file parameters.
	curr_datetime = datetime('now');

	time_of_day_str = char(timeofday( curr_datetime ));
	colon_locs = find(time_of_day_str == ':');
	time_of_day_str(colon_locs) = '_';

	date_str = [ num2str( month( curr_datetime ) ) '_' num2str( day( curr_datetime ) ) ];

	if ~exist('save_file_name')
		save_file_name = ['data/ecc/oc74_interm_results_' num2str(c_sq.dim_x) 'x' num2str(c_sq.dim_y) 'drones_' ...
							date_str '_' time_of_day_str '.mat'];
	end

	run_opt_flag = true;

	if load_data_flag
		load(save_file_name)
		load_data_flag = true;
	end

	ops = sdpsettings('verbose',verbosity);

	results_filename_base = 'results/ecc2020/oc74_results_';

	drone_img_side = 0.6;

	%%%%%%%%%%%
	%% Tests %%
	%%%%%%%%%%%

	target_set_factor = 0.5;
	target_set = Px0*target_set_factor;

	if ~exist('BG')
		[BG,contr,opt_data,BG_timing] = in_sys.synth_robust_reach_contr(Px0,Pu,'P_target',target_set, 'debug' , verbosity)
		save(save_file_name,'contr','opt_data','BG','c_sq','eta_u','eta_x','BG_timing','target_set','Px0','curr_datetime')
	elseif ~exist('contr')
		[~,contr,opt_data] = in_sys.synth_robust_reach_contr(Px0,Pu,'P_target',target_set,'BG',BG, 'debug' , verbosity)
		save(save_file_name,'contr','opt_data','BG','c_sq','eta_u','eta_x','BG_timing','target_set','Px0')
	else
		warning('This indicates that you are only interested in plotting.')
	end

	results.BG = BG;
	results.BGAlgorithmTime = BG_timing;
	results.contr = contr;
	results.opt_data = opt_data;

	figure;
	BG.plot();
	saveas(gcf,[results_filename_base 'belief_graph'], 'epsc')


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Simulate the System With this Form of Controller %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%rand_word_idx = randi( BG.ModeLanguage.cardinality() );
	worst_word_idx = 1;

	sig = BG.ModeLanguage.words{worst_word_idx};

	in_w = unifrnd(-eta_w,eta_w,2,length(sig));
	for symb_idx = 1:length(sig)
		if sig(symb_idx) == 2
			in_w(:,symb_idx) = 1.2*eta_w*ones(2,1);
		end
	end

	in_x0 = [-eta_x;eta_x;eta_x;-eta_x];

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

	M2_prime = max([M2 + 0.3,max_x+0.3,max_y+0.3]);
	axis_lims = [ 	max(-M2_prime,min_x-0.3) , M2_prime , ...
					max(-M2_prime,min_x-0.3) , M2_prime ];


	if length(sig) < 4
		error('This plot was designed for languages of length at least equal to 4.')
	end

	figure;
	hold on;
	%set(gca,'color','none')
	plot_idx = 0;
	if length(sig) == 4
		tau_list = [0,1,2,3,4];
	elseif length(sig) > 4
		tau_list = [0,1,2, length(sig)-2 ,length(sig)-1 ,length(sig)];
	end
		
	for tau = tau_list
		plot_idx = plot_idx + 1;
		subplot(1,length(tau_list),plot_idx)
		hold on;

		time_idx = tau+1;

		%Plot Bounding Set (Initial State Set)
		rectangle('Position',[center_pos(1)-M2,center_pos(2)-M2, ...
										2*M2, ...
										2*M2 ],...
							'FaceColor','none');
		%Plot Target Set
		rectangle('Position',[center_pos(1)-target_set_factor*eta_x,center_pos(2)-target_set_factor*eta_x, ...
										2*target_set_factor*eta_x, ...
										2*target_set_factor*eta_x ],...
							'FaceColor',[0.3010 0.7450 0.9330],'LineStyle', '--')

		%Plot Robots
		dx = drone_img_side ; dy = drone_img_side ;              %# Add to the plot
		x1min = x(1,time_idx)-dx ; x1max = x(1,time_idx)+dx ;
		y1min = x(2,time_idx)-dy ; y1max = x(2,time_idx)+dy ;
		x2min = x(3,time_idx)-dx ; x2max = x(3,time_idx)+dx ;
		y2min = x(4,time_idx)-dy ; y2max = x(4,time_idx)+dy ;
		% Make background transperent
		[img,map,alphachannel] = imread('data/image88.png');
		%img = flipud(img) ;

		disp(['Are any of the pixels in image88 not 0 or 1?' num2str(any(any( (img ~= 0) & (img~=1) ))) ])
		size(img)

		h1 = image([x1min x1max],[y1min y1max],img,'AlphaData',alphachannel);        % Plot the robot image
		h2 = image([x2min x2max],[y2min y2max],img,'AlphaData',alphachannel);        % Plot the robot image
		
		% scatter(x(1,time_idx),x(2,time_idx)) %Plot Robot 1
		% scatter(x(3,time_idx),x(4,time_idx)) %Plot Robot 2

		axis(axis_lims)
		title(['$t=' num2str(tau) '$'],'Interpreter','latex')
	end
	set(gcf,'units','Normalized','Position',[0 0 1.0 0.7])

	saveas(gcf,[results_filename_base 'comic_strip'], 'epsc')

	%% Plot Individual Frames %%

	plot_idx = 0;
	if length(sig) == 4
		tau_list = [0,1,2,3,4];
	elseif length(sig) > 4
		tau_list = [0,1,2, length(sig)-2 ,length(sig)-1 ,length(sig)];
	end
		
	% for tau = tau_list
	% 	figure;
	% 	hold on;
	% 	set(gca,'color','none')

	% 	plot_idx = plot_idx + 1;

	% 	time_idx = tau+1;

	% 	%Plot Bounding Set (Initial State Set)
	% 	rectangle('Position',[center_pos(1)-M2,center_pos(2)-M2, ...
	% 									2*M2, ...
	% 									2*M2 ],...
	% 						'FaceColor','none');
	% 	%Plot Target Set
	% 	rectangle('Position',[center_pos(1)-target_set_factor*eta_x,center_pos(2)-target_set_factor*eta_x, ...
	% 									2*target_set_factor*eta_x, ...
	% 									2*target_set_factor*eta_x ],...
	% 						'FaceColor',[0.3010 0.7450 0.9330],'LineStyle', '--')

	% 	%Plot Robots
	% 	dx = drone_img_side ; dy = drone_img_side ;              %# Add to the plot
	% 	x1min = x(1,time_idx)-dx ; x1max = x(1,time_idx)+dx ;
	% 	y1min = x(2,time_idx)-dy ; y1max = x(2,time_idx)+dy ;
	% 	x2min = x(3,time_idx)-dx ; x2max = x(3,time_idx)+dx ;
	% 	y2min = x(4,time_idx)-dy ; y2max = x(4,time_idx)+dy ;
	% 	% Make background transperent
	% 	[img,map,alphachannel] = imread('data/image88.png');
	% 	%img = flipud(img) ;

	% 	disp(['Are any of the pixels in image88 not 0 or 1?' num2str(any(any( (img ~= 0) & (img~=1) ))) ])
	% 	size(img)

	% 	h1 = image([x1min x1max],[y1min y1max],img,'AlphaData',alphachannel);        % Plot the robot image
	% 	h2 = image([x2min x2max],[y2min y2max],img,'AlphaData',alphachannel);        % Plot the robot image
		
	% 	% scatter(x(1,time_idx),x(2,time_idx)) %Plot Robot 1
	% 	% scatter(x(3,time_idx),x(4,time_idx)) %Plot Robot 2

	% 	axis([-M2_prime, M2_prime, -M2_prime, M2_prime])
	% 	title(['$t=' num2str(tau) '$'],'Interpreter','latex')

	% 	set(gcf,'units','Normalized','Position',[0 0 0.9 0.6])

	% 	saveas(gcf,[results_filename_base '_frame' num2str(tau) ], 'epsc')
	% end
	

	results.x = x;
	results.u = u;
	results.y = y;
	results.sig = sig;
	results.b = contr.b_hist

	results.temp_img = img;

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