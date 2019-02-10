function [results] = observer_comparison37(varargin)
	%observer_comparison37.m
	%	Verifying the behavior of a simple Lane Keeping system.

	%%%%%%%%%%%%%%%
	%% Constnats %%
	%%%%%%%%%%%%%%%

	dt = 0.01;

	%Missing Data Parameters
	word_len = 6;
	L = {};
	for L_ind = 2:word_len-1
		L{L_ind-1} = ones(1,word_len);
		L{L_ind}(L_ind) = 0;
	end
	L{word_len-1} = ones(1,word_len);

	%Simulation constants
	num_runs = 100;

	%target sets
	M1 = 0.3;
	M2 = 0.5;

	%System Parameters
	A = [ zeros(2,1) [1;-20] ];
	n_x = size(A,1);

	B = [ 0; 1];
	F = zeros(n_x,1);
	
	C = eye(2);%[1,0];
	n_y = size(C,1);

	temp_sys = ss(A,B,C,0);
	dt = 0.05;
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

	%Simulate
	[ ctrl_sim1, ~ ] = oc37_contr1.simulate_n_runs( simple_LK , M1 , num_runs )

	ctrl_sim1_mod = [];
	for ind = 1:length(ctrl_sim1)
		ctrl_sim1_mod(:,:,ind) = simple_LK.C*ctrl_sim1{ind};
	end

	figure;
	hold on;
	bar([0:word_len]*dt,[M1 ones(1,word_len-1)*M2 M1],1,'w')
	bar([0:word_len]*dt,-1*[M1 ones(1,word_len-1)*M2 M1],1,'w')
	plot([0:word_len]*dt,(lane_width/2)*ones(1,word_len+1),'k','LineWidth',2)
	plot([0:word_len]*dt,-(lane_width/2)*ones(1,word_len+1),'k','LineWidth',2)
	plot([0:word_len]*dt,zeros(1,word_len+1),'k--','LineWidth',2)


	for sim_num = 1:10
		plot([0:word_len]*dt,ctrl_sim1_mod(1,:,sim_num))
	end

	axis([0,word_len*dt,-1,1])

	% Saving

	results.exp1.sys = simple_LK;
	results.exp1.opt = oc37_opt1;
	results.exp1.contr = oc37_contr1;
end