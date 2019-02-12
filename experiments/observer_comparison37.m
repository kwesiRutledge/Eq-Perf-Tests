function [results] = observer_comparison37(varargin)
	%observer_comparison37.m
	%	Verifying the behavior of a simple Lane Keeping system.

	%%%%%%%%%%%%%%%
	%% Constnats %%
	%%%%%%%%%%%%%%%

	%Missing Data Parameters
	word_len = 12;
	md_window = [1 10];
	L = {};
	for L_ind = md_window(1):md_window(2)
		L{end+1} = ones(1,word_len);
		L{end}(L_ind) = 0;
		L{end}(L_ind+1) = 0;
	end
	L{end+1} = ones(1,word_len);

	% L = {ones(1,word_len)};

	%Simulation constants
	num_runs = 100;

	%target sets
	M1 = 0.3;
	% M2 = 0.5;

	%System Parameters
	A = [ zeros(2,1) [1;-20] ];
	n_x = size(A,1);

	B = [ 0; 1];
	F = zeros(n_x,1);
	
	C = eye(2); %[1 0];
	n_y = size(C,1);

	temp_sys = ss(A,B,C,0);
	dt = 0.1;
	temp_dsys = c2d(temp_sys,dt);

	temp_sys2 = ss(A,[1;0],C,0);
	temp_dsys2 = c2d(temp_sys2,dt);

	eta_w = 0.05; eta_v = 0.1;

	simple_LK = Aff_Dyn(	temp_dsys.A,temp_dsys.B,F,C,...
							eta_w,eta_v, ...
							temp_dsys2.B, eye(n_y) );

	lane_width = 2*0.9;

	%%%%%%%%%%%%%%%
	%% Synthesis %%
	%%%%%%%%%%%%%%%

	[ oc37_opt1 , oc37_contr1 ] = simple_LK.eq_rec_design_pb( 'Min_M2' , M1 , L );

	%%%%%%%%%%%%%
	%% Results %%
	%%%%%%%%%%%%%

	M2 = oc37_opt1.M2;

	runs_per_word = 1;
	fs = 20;

	%Simulate
	% ctrl_sim1 = {};
	% for L_ind = 1: length(L)
	% 	[ ctrl_sim1{end+1}, ~ ] = oc37_contr1.simulate_n_runs( simple_LK , M1 , runs_per_word , 'in_sigma' , L{L_ind} , 'in_w' , eta_w*ones(1,word_len) );
	% end

	ad = simple_LK;
	contr = oc37_contr1;

	x0_og = [M1;M1];
	w_og = eta_w*ones(1,word_len);
	v_og = [0.0038,  0.0657,  -0.0013, -0.0249, -0.0030, -0.0056, 0.0550,  0.0072, -0.0829, 0.0078, -0.0264, 0.0265;
	        -0.0841, -0.0454, -0.0797, -0.0543, -0.0031, 0.0907,  -0.0695, 0.0446, -0.0759, 0.0410, -0.0807, 0.0022];
	v_og2 = -eta_v*ones(2,word_len);
	%Constants
	for L_ind = 1:length(L)
		sig = L{L_ind};
		T = length(sig);
		for run_ind = 1:1
			% T = size(obj.L,2);
			n = size(ad.A,1);
			wd = size(ad.B_w,2);
			vd = size(ad.C_v,2);

			%Generate Random Variables
			x0 = x0_og;
			w  = w_og;
			v  = v_og2;

			%Simulate system forward.
			x_t = x0;
			y_t = ad.C*x0 + ad.C_v*v(:,1);
			
			x_0_t = x_t;
			y_0_t = y_t;
			x_tp1 = Inf(n,1);
			sig = [sig,1];
			for t = 0:T-1
				%Use Affine Dynamics with proper control law.
				x_tp1 = ad.A * x_t + ...
						ad.B * contr.apply_control( sig([1:t+1]) , y_0_t ) + ...
						ad.B_w * w(:,t+1) + ...
						ad.f ;
				%Update other variables in system
				x_t = x_tp1;
				x_0_t = [x_0_t x_t];

				if( t == T-1 )
					continue;
				end

				y_t = sig(t+2)*(ad.C*x_t + ad.C_v*v(:,t+1));
				y_0_t = [y_0_t; y_t];
			end
		end
		ctrl_sim1{L_ind} = {x_0_t};
	end

	ctrl_sim1_mod = [];
	for ind = 1:length(ctrl_sim1)
		for ind2 = 1:length(ctrl_sim1{ind})
			ctrl_sim1_mod(:,:,(ind-1)*length(ctrl_sim1{ind})+ind2) = simple_LK.C*ctrl_sim1{ind}{ind2};
		end
	end

	figure;
	hold on;

	%
	road_center_line = kron([1:word_len+1]-1,[dt;dt*eta_w]);

	% bar([0:word_len]*dt,[M1 ones(1,word_len-1)*M2 M1],1,'w')
	% bar([0:word_len]*dt,-1*[M1 ones(1,word_len-1)*M2 M1],1,'w')
	%Create series of clear boxes for M1 guarantees
	for t = 1:word_len+1
		if (t == 1) || (t==word_len+1)
			rect_dim = [dt 2*M1];
		else
			rect_dim = [dt 2*M2];
		end
		%Adjust rectangle position so that it is centered where we want it.
		rect_xy = road_center_line(:,t) - (1/2)*rect_dim';
		rectangle('Position',[rect_xy' rect_dim])
	end
	plot(road_center_line(1,:),road_center_line(2,:)+(lane_width/2)*ones(1,word_len+1),'k','LineWidth',2)
	plot(road_center_line(1,:),road_center_line(2,:)-(lane_width/2)*ones(1,word_len+1),'k','LineWidth',2)
	plot(road_center_line(1,:),road_center_line(2,:),'k--','LineWidth',2)

	for sim_num = 1:size(ctrl_sim1_mod,3)
		plot([0:word_len]*dt,ctrl_sim1_mod(1,:,sim_num))
	end

	axis([0,word_len*dt,min(road_center_line(2,:)-(lane_width/2)*ones(1,word_len+1)),max(road_center_line(2,:)+(lane_width/2)*ones(1,word_len+1))])
	xlabel('Time $t$ (s)','Interpreter','latex','Fontsize',fs)
	ylabel('Deviation from Lane Center, $x_1(t)$','Interpreter','latex','Fontsize',fs)
	legend('Road Boundary','Road Boundary','Road Centerline')

	% Saving

	results.exp1.sys = simple_LK;
	results.exp1.M1  = M1;
	results.exp1.L   = L;
	results.exp1.opt = oc37_opt1;
	results.exp1.contr = oc37_contr1;
	results.exp1.trajs = ctrl_sim1_mod;
end