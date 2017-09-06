% acc_finite_horizon_control.m

clear all;
close all;
clc;

%% Add Functions to Path
addpath('./functions')

%% Load ACC System
% load('../data/acc_cont_n_discr_params.mat')

% % Create Structs that hold information about Continuous and Discrete Systems
% acc_dsys.A = A_d; acc_dsys.B = B_d; acc_dsys.C = C_d; acc_dsys.E = E_d; acc_dsys.F = F_d;
% acc_csys.A = A_c; acc_csys.B = B_c; acc_csys.C = C_c; acc_csys.E = E_c; acc_csys.F = F_c;

% % Save these structs for future use
% save('../data/acc_m1_systems.mat','acc_csys','acc_dsys');

load('../data/acc_m1_systems.mat')

%% Constants

T = 50;

eps = 10e-4;

perform_experiment = [
false;	%1
false;	%2
false;	%3
false;	%4
false;	%5
false;	%6
false;	%7
true;	%8
false;	%9
false;	%10
true 	%11
];

%% Experiment 1
disp('============================================')
disp('Experiment 1: Basics of Skaf and Boyd Method')

if perform_experiment(1)

	acc_dsys.x0 = [0;1;0];
	acc_dsys.d = 0.1;
	acc_dsys.m = 0.05;
	n = size(acc_dsys.A,1);

	% Create Matrices Necessary for Original Skaf and Boyd Method
	%------------------------------------------------------------
	[G_e1,H_e1,Cm_e1,x0m_e1] = create_skaf_n_boyd_matrices(acc_dsys,T);

	disp('Created Skaf and Boyd Matrices.')

	w_e1 = unifrnd(-acc_dsys.d,acc_dsys.d,n,T);
	v_e1 = unifrnd(-acc_dsys.m,acc_dsys.m,size(acc_dsys.C,1),T);

	%Reshape w and v into columns
	w_e1_rs	= reshape(w_e1,prod(size(w_e1)),1);
	v_e1_rs = reshape(v_e1,prod(size(v_e1)),1); 

	cvx_begin quiet
		variable Q(size(H_e1,2),size(Cm_e1,1))
		variable r(size(H_e1,2))

		Pxw = (eye(n*(T+1))+H_e1*Q*Cm_e1)*G_e1;
		Pxv = H_e1*Q;
		x_tilde = (eye(n*(T+1)) + H_e1*Q*Cm_e1)*x0m_e1 + H_e1*r;

		R = [zeros(n,n*T) eye(n)];

		minimize norm( R*( x_tilde + Pxw * w_e1_rs + Pxv * v_e1_rs ),Inf)
		subject to
			% Q (n,n) block lower triangular
		    for i=0:T-1 
		        Q(i+1,i+2:end) == 0
		    end
	cvx_end

	%Plot the results of the found controller
	e1_results.x = reshape( x_tilde + Pxw * w_e1_rs + Pxv * v_e1_rs , n , length(x_tilde)/n );
	for t = 1:size(e1_results.x,2)
		e1_results.x_inf_norm(t) = norm(e1_results.x(:,t),Inf);
	end

	disp('Displaying how the designed controller''s infinity norm is handled throughout the time window.')

	figure(1);
	plot(e1_results.x_inf_norm)
	title(['Infinity Norm of X in T=' num2str(T) ' Window'])
	xlabel('Time Index')
	ylabel('$||x||_{\infty}$','Interpreter','latex')
	axis([0 60 0 max(e1_results.x_inf_norm)+0.5])

	%How well do these gains generalize?
	num_tests = 300;
	w2_e1 = unifrnd(-acc_dsys.d,acc_dsys.d,n,T,num_tests);
	v2_e1 = unifrnd(-acc_dsys.m,acc_dsys.m,size(acc_dsys.C,1),T,num_tests);

	x_test = [];

	disp('Testing how well the randomly generated noises from 1 trajectory generalize')
	disp(' to how well the infinity norm is controlled in general.')

	%Test each one using the same feedback matrices.
	inf_norm_traj = [];
	for test_num = 1 : num_tests
		%Create temporary w and v vectors.
		temp_w = w2_e1(:,:,test_num);
		temp_v = v2_e1(:,:,test_num);

		temp_w_rs = reshape(temp_w,prod(size(temp_w)),1);
		temp_v_rs = reshape(temp_v,prod(size(temp_v)),1);

		%Apply same feedback matrix/gains to calculate x.
		temp_x = x_tilde + Pxw * temp_w_rs + Pxv * temp_v_rs;

		%Calculate Infinity Norm Trajectory
		temp_x_rs = reshape(temp_x,n,length(temp_x)/n);

		for col = 1:size(temp_x_rs,2)
			inf_norm_traj(test_num,col) = norm(temp_x_rs(:,col),Inf);
		end

	end

	avg_inf_norm_traj = mean(inf_norm_traj);

	figure(2);
	plot(avg_inf_norm_traj)
	title('Average Trajectory of $||x||_{\infty}$ on Test Data','Interpreter','latex')
	xlabel('Time Index')
	ylabel('$||x||_{\infty}$','Interpreter','latex')
	axis([0 60 0 max(e1_results.x_inf_norm)+0.5])

else
	disp('User opted out of Experiment 1.')
end

%% Experiment 2

disp('============================================================')
disp('Experiment 2: Comparing Skaf and Boyd against Yong''s Method')
disp('For the same dataset when minimizing the norm error at the Tth timestep')

if perform_experiment(2)

	acc_dsys.x0 = [0;1;0];
	acc_dsys.d = 0.1;
	acc_dsys.m = 0.05;
	n = size(acc_dsys.A,1);

	% Create Matrices Necessary for Original Skaf and Boyd Method
	%------------------------------------------------------------
	[G_e2_skaf,H_e2_skaf,Cm_e2_skaf,x0m_e2_skaf] = create_skaf_n_boyd_matrices(acc_dsys,T);

	disp('Created Skaf and Boyd Matrices.')

	% Create many different disturbances to train/create our "Q" matrix.
	%-----------------------------------------------------------------------
	num_train_examples_e2 = 5;

	w_train_e2 = unifrnd(-acc_dsys.d,acc_dsys.d,n,T,num_train_examples_e2);
	v_train_e2 = unifrnd(-acc_dsys.m,acc_dsys.m,size(acc_dsys.C,1),T,num_train_examples_e2);

	disp('Created Disturbances for training.')

	% Create Q and r Matrices
	%------------------------
	w_train_e2_rs = reshape(w_train_e2,size(w_train_e2,1)*size(w_train_e2,2),1,size(w_train_e2,3));
	v_train_e2_rs = reshape(v_train_e2,size(v_train_e2,1)*size(v_train_e2,2),1,size(v_train_e2,3));

	cvx_begin quiet
		variable Q(size(H_e2_skaf,2),size(Cm_e2_skaf,1))
		variable r(size(H_e2_skaf,2))

		Pxw = (eye(n*(T+1))+H_e2_skaf*Q*Cm_e2_skaf)*G_e2_skaf;
		Pxv = H_e2_skaf*Q;
		x_tilde = (eye(n*(T+1)) + H_e2_skaf*Q*Cm_e2_skaf)*x0m_e2_skaf + H_e2_skaf*r;

		R = [zeros(n,n*T) eye(n)];

		%Create vector containing ALL of the infinity norms of our examples
		for i = 1 : num_train_examples_e2
			y(i) = norm( R*( x_tilde + Pxw * w_train_e2_rs(:,:,i) + Pxv * v_train_e2_rs(:,:,i) ),Inf);
		end

		minimize y*ones(size(y'))*(1/num_train_examples_e2)
		subject to
			% Q (n,n) block lower triangular
		    for i=0:T-1 
		        Q(i+1,i+2:end) == 0
		    end
	cvx_end

	if strcmp(cvx_status,'Solved')
		disp('Q and r matrices found through optimization')
	else
		error('Optimization to find Q and r failed.')
	end

	% Plot the results of the found controller
	%-----------------------------------------

	e2_results.skaf_training_optval = cvx_optval;
	e2_results.skaf_Q = Q;
	e2_results.skaf_r = r;

	avg_norm_traj_e2_skaf = zeros(1,T+1);
	for example_num = 1 : num_train_examples_e2
		
		%Reconstruct X using our Q
		temp_x = x_tilde + Pxw * w_train_e2_rs(:,:,example_num) + Pxv * v_train_e2_rs(:,:,example_num);
		for time_ind = 1 : T+1
			temp_x_inf_norm(time_ind) = norm( temp_x([(time_ind-1)*n+1:time_ind*n]) ,Inf);
		end
		avg_norm_traj_e2_skaf = avg_norm_traj_e2_skaf + (1/num_train_examples_e2)*temp_x_inf_norm;
	end

	e2_results.skaf_avg_norm_traj = avg_norm_traj_e2_skaf;

	figure;
	plot(avg_norm_traj_e2_skaf)
	title('E2: Boyd and Skaf Method Used to Minimize $||x(T)||_{\infty}$','Interpreter','latex')
	xlabel('Time Index t')
	ylabel('Average $||x(t)||_{\infty}$','Interpreter','latex')

	disp('Skaf and Boyd Performance Noted.')

	% Create Yong Modified Matrices
	%------------------------------
	s = n; %Size of the output state is the same as the actual state

	acc_dsys_dyn_out.A = [ acc_dsys.A zeros(n,s) ; zeros(s,n) zeros(s) ];
	acc_dsys_dyn_out.B = [ zeros(n,s) acc_dsys.B ; eye(s) zeros(s,size(acc_dsys.B,2)) ];
	acc_dsys_dyn_out.C = [ zeros(s,n) eye(s) ; acc_dsys.C zeros(size(acc_dsys.C,1),s) ];

	acc_dsys_dyn_out.x0 = [ acc_dsys.x0 ; zeros(size(acc_dsys.x0)) ];
	acc_dsys_dyn_out.d = 0.1;
	acc_dsys_dyn_out.m = 0.05;

	[G_e2_yong,H_e2_yong,Cm_e2_yong,x0m_e2_yong] = create_skaf_n_boyd_matrices(acc_dsys_dyn_out,T);

	%Recreate w and v matrices
	% w_bar_train_e2 = []; v_bar_train_e2 = [];
	for example_num = 1 : num_train_examples_e2
			% w_bar_train_e2(:,:,example_num) = [];
			% v_bar_train_e2(:,:,example_num) = [];
			for time_ind = 1:T
				w_bar_train_e2(:,time_ind,example_num) = [ eye(size(w_train_e2,1)) ; zeros(size(w_train_e2,1)) ]*w_train_e2(:,time_ind);
				v_bar_train_e2(:,time_ind,example_num) = [ zeros(s,size(acc_dsys.C,1)) ; eye(size(acc_dsys.C,1)) ]*v_train_e2(:,time_ind);
			end
	end

	disp('Created Yong Matrices.')

	% Create F_star and r_star
	%-------------------------
	w_bar_train_e2_rs = reshape(w_bar_train_e2,size(w_bar_train_e2,1)*size(w_bar_train_e2,2),1,size(w_bar_train_e2,3));
	v_bar_train_e2_rs = reshape(v_bar_train_e2,size(v_bar_train_e2,1)*size(v_bar_train_e2,2),1,size(v_bar_train_e2,3));

	cvx_begin quiet
		variable Q_star(size(H_e2_yong,2),size(Cm_e2_yong,1))
		variable r_star(size(H_e2_yong,2))

		Pxw = (eye(2*n*(T+1))+H_e2_yong*Q_star*Cm_e2_yong)*G_e2_yong;
		Pxv = H_e2_yong*Q_star;
		x_tilde = (eye(2*n*(T+1)) + H_e2_yong*Q_star*Cm_e2_yong)*x0m_e2_yong + H_e2_yong*r_star;

		R = [zeros(n,(n+s)*T) eye(n) zeros(n,s)];

		%Create vector containing ALL of the infinity norms of our examples
		for i = 1 : num_train_examples_e2
			y2(i) = norm( R*( x_tilde + Pxw * w_bar_train_e2_rs(:,:,i) + Pxv * v_bar_train_e2_rs(:,:,i) ),Inf);
		end

		minimize y2*ones(size(y2'))*(1/num_train_examples_e2)
		subject to
			% Q_star (2n,2n) block lower triangular
		    for i=0:T-1 
		        Q_star([2*i+1:2*i+2],2*i+3:end) == 0
		    end
	cvx_end

	if strcmp(cvx_status,'Solved')
		disp('Created Q_star and r_star through optimization.')
	else
		disp('Optimization for Q_star and r_star failed.')
	end

	%Plot the Results of Yong's Dynamic Output Feedback Controller
	e2_results.yong_training_optval = cvx_optval;
	e2_results.yong_Q = Q_star;
	e2_results.yong_r = r_star;

	avg_norm_traj_e2_yong = zeros(1,T+1);
	for example_num = 1 : num_train_examples_e2
		
		%Reconstruct X using our Q
		temp_x = x_tilde + Pxw * w_bar_train_e2_rs(:,:,example_num) + Pxv * v_bar_train_e2_rs(:,:,example_num);
		for time_ind = 1 : T+1
			temp_x_inf_norm(time_ind) = norm( temp_x([(time_ind-1)*(n+s)+1:time_ind*(n+s) - s]) ,Inf);
		end
		avg_norm_traj_e2_yong = avg_norm_traj_e2_yong + (1/num_train_examples_e2)*temp_x_inf_norm;
	end

	e2_results.yong_avg_norm_traj = avg_norm_traj_e2_yong;

	figure;
	plot(avg_norm_traj_e2_yong)
	title('E2: Yong Method Used to Minimize $||x(T)||_{\infty}$','Interpreter','latex')
	xlabel('Time Index t')
	ylabel('Average $||x(t)||_{\infty}$','Interpreter','latex')

	figure;
	hold on;
	plot(avg_norm_traj_e2_yong)
	plot(avg_norm_traj_e2_skaf,':')
	legend('Yong','Skaf')
	xlabel('Time Index t')
	ylabel('Average $||x(t)||_{\infty}$','Interpreter','latex')
	title('E2: Yong and Skaf Methods compared, when minimizing $||x(T)||_{\infty}$','Interpreter','latex')

	% Save Data to Data folder
	save('data/acc_skaf_v_yong_comparison_e2.mat','e2_results')

else
	disp('Experiment not performed.')
end

%% Experiment 3
if perform_experiment(3)

	%Announce Experiment
	disp('================================================')
	disp('Experiment 3: Comparison of Robust Optimizations')
	disp(' ');


	%Load data
	if exist('acc_dsys')
		disp('Data already loaded.')
	else
		load('../data/acc_m1_systems.mat')
	end

	%Defining Constants
	acc_dsys.x0 = [0;1;0];
	acc_dsys.d = 0.1;
	acc_dsys.m = 0.05;
	n = size(acc_dsys.A,1);

	e3_params.n = n;
	e3_params.d = acc_dsys.d;
	e3_params.m = acc_dsys.m;
	e3_params.T = 3;

	% Create Matrices Necessary for Original Skaf and Boyd Method
	%------------------------------------------------------------
	[G_e3,H_e3,Cm_e3,x0m_e3] = create_skaf_n_boyd_matrices(acc_dsys,e3_params.T);
	disp('Created Skaf and Boyd Matrices.')

	% Perform Robust Optimization Using YALMIP
	%-----------------------------------------
	e3_params.Q_yalmip = sdpvar(size(H_e3,2),size(Cm_e3,1),'full');
	e3_params.r_yalmip = sdpvar(size(H_e3,2),1,'full');

	%Disturbance vectors
	w_e3_yalmip = sdpvar(n*e3_params.T,1,'full');
	v_e3_yalmip = sdpvar(size(acc_dsys.C,1)*e3_params.T,1,'full');
	disp('Created YALMIP Optimization Variables.')

	%Create Expressions Containing Q,r for Optimization
	Pxw = (eye(n*(e3_params.T+1))+H_e3*e3_params.Q_yalmip*Cm_e3)*G_e3;
	Pxv = H_e3*e3_params.Q_yalmip;
	x_tilde = (eye(n*(e3_params.T+1)) + H_e3*e3_params.Q_yalmip*Cm_e3)*x0m_e3 + H_e3*e3_params.r_yalmip;

	%Create Objective
	R = [zeros(n,n*e3_params.T) eye(n)]; %Create selection matrix
	e3_params.obj_yalmip = norm( R*(x_tilde + Pxw * w_e3_yalmip + Pxv * v_e3_yalmip) , Inf );
	disp('Created YALMIP Objective.')

	%Create Constraints

	%Q is lower diagonal.
	l_diag_constr = [];
	for bl_row_num = 1 : e3_params.T-1
		l_diag_constr = l_diag_constr + [ e3_params.Q_yalmip(	[(bl_row_num-1)*size(acc_dsys.B,2)+1:bl_row_num*size(acc_dsys.B,2)], ...
																[bl_row_num*size(acc_dsys.C,1)+1:end] ) == 0 ];
    end

    l_diag_constr

    %Robustifying against w and v
    robust_constr = [];
    robust_constr = robust_constr + [ -e3_params.m <= v_e3_yalmip <= e3_params.m , uncertain(v_e3_yalmip) ];
    robust_constr = robust_constr + [ -e3_params.d <= w_e3_yalmip <= e3_params.d , uncertain(w_e3_yalmip) ]

    disp('Created YALMIP Constraints.')

    %Solve Optimization
    disp('Solving YALMIP Robust Optimization...')
    ops = sdpsettings('verbose',0);
    e3_results.sol_yalmip = optimize(l_diag_constr+robust_constr,e3_params.obj_yalmip,ops);

    if e3_results.sol_yalmip.problem == 0
    	disp('YALMIP Optimization Solved')
    else
    	error('YALMIP Optimization NOT Solved.')
    end

    e3_results.Q_yalmip = value(e3_params.Q_yalmip);
    e3_results.r_yalmip = value(e3_params.r_yalmip);

    clear Pxw Pxv x_tilde

    % Perform Robust Optimization Using Duality and CVX
    %--------------------------------------------------
    disp(' ')

    error('The CVX Optimization does not appear to be solveable.')

    cvx_begin
    	variable Q_e3_cvx(size(H_e3,2),size(Cm_e3,1))
    	variable r_e3_cvx(size(H_e3,2))

    	Pxw = (eye(n*(e3_params.T+1))+H_e3*Q_e3_cvx*Cm_e3)*G_e3;
		Pxv = H_e3*Q_e3_cvx;
		x_tilde = (eye(n*(e3_params.T+1)) + H_e3*Q_e3_cvx*Cm_e3)*x0m_e3 + H_e3*r_e3_cvx;

		variables alpha_e3(1) lambda_e3(size(R,1),2) eta_e3(size(Pxw,2),2) gamma_e3(size(R,1),2) theta_e3(size(Pxv,2),2)

		disp('Created CVX Variables.')

		minimize alpha_e3
		subject to
			%Lower Diagonal Constraint on Q
			for bl_row_num = 1 : e3_params.T-1
				Q_e3_cvx(	[(bl_row_num-1)*size(acc_dsys.B,2)+1:bl_row_num*size(acc_dsys.B,2)], ...
							[bl_row_num*size(acc_dsys.C,1)+1:end] ) == 0
    		end
    		%Constraints coming from transformation of Robust Optimization problem
    		-(lambda_e3(:,1) - lambda_e3(:,2))'*R*x_tilde ...
    			- e3_params.d*sum( eta_e3(:,1) + eta_e3(:,2) ) ...
    			- e3_params.m * sum( theta_e3(:,1) + theta_e3(:,2) ) <= alpha_e3
    		sum(gamma_e3(:,1) + gamma_e3(:,2)) == 1
    		lambda_e3(:,1) + lambda_e3(:,2) == -gamma_e3(:,1)+gamma_e3(:,2)
    		(-lambda_e3(:,1) + lambda_e3(:,2))'*R*Pxw == (eta_e3(:,1) - eta_e3(:,2))'
    		(-lambda_e3(:,1) + lambda_e3(:,2))'*R*Pxv == (theta_e3(:,1) - theta_e3(:,2))'
    		lambda_e3 >= 0
    		theta_e3 >= 0
    		gamma_e3 >= 0
    		theta_e3 >= 0
    cvx_end

    e3_results.Q_cvx = Q_e3_cvx;
    e3_results.r_cvx = r_e3_cvx;


else
	disp('Experiment not performed.')
end

%% Experiment 4: How low can the final error become when T varies

if perform_experiment(4)

	%Announce Experiment
	disp('================================================')
	disp('Experiment 4: Robust Optimizations, Varying T')
	disp(' ');


	%Load data
	if exist('acc_dsys')
		disp('Data already loaded.')
	else
		load('../data/acc_m1_systems.mat')
	end

	%Defining Constants
	acc_dsys.x0 = [0;1;0];
	acc_dsys.d = 0.1;
	acc_dsys.m = 0.05;
	n = size(acc_dsys.A,1);

	e4_params.n = n;
	e4_params.d = acc_dsys.d;
	e4_params.m = acc_dsys.m;
	e4_params.T_max = 20;

	for T_e4 = 1 : e4_params.T_max


		% Create Matrices Necessary for Original Skaf and Boyd Method
		%------------------------------------------------------------
		[G_e4,H_e4,Cm_e4,x0m_e4] = create_skaf_n_boyd_matrices(acc_dsys,T_e4);
		disp('Created Skaf and Boyd Matrices.')

		% Perform Robust Optimization Using YALMIP
		%-----------------------------------------
		e4_params.Q = sdpvar(size(H_e4,2),size(Cm_e4,1),'full');
		e4_params.r = sdpvar(size(H_e4,2),1,'full');

		%Disturbance vectors
		w_e4_yalmip = sdpvar(n*T_e4,1,'full');
		v_e4_yalmip = sdpvar(size(acc_dsys.C,1)*T_e4,1,'full');
		disp('Created YALMIP Optimization Variables.')

		%Create Expressions Containing Q,r for Optimization
		Pxw = (eye(n*(T_e4+1))+H_e4*e4_params.Q*Cm_e4)*G_e4;
		Pxv = H_e4*e4_params.Q;
		x_tilde = (eye(n*(T_e4+1)) + H_e4*e4_params.Q*Cm_e4)*x0m_e4 + H_e4*e4_params.r;

		%Create Objective
		R = [zeros(n,n*T_e4) eye(n)]; %Create selection matrix
		e4_params.obj_yalmip = norm( R*(x_tilde + Pxw * w_e4_yalmip + Pxv * v_e4_yalmip) , Inf );
		disp('Created YALMIP Objective.')

		%Create Constraints

		%Q is lower diagonal.
		l_diag_constr = [];
		for bl_row_num = 1 : e3_params.T-1
			l_diag_constr = l_diag_constr + [ e3_params.Q_yalmip(	[(bl_row_num-1)*size(acc_dsys.B,2)+1:bl_row_num*size(acc_dsys.B,2)], ...
																	[bl_row_num*size(acc_dsys.C,1)+1:end] ) == 0 ];
	    end

	    l_diag_constr

	    %Robustifying against w and v
	    robust_constr = [];
	    robust_constr = robust_constr + [ -e3_params.m <= v_e3_yalmip <= e3_params.m , uncertain(v_e3_yalmip) ];
	    robust_constr = robust_constr + [ -e3_params.d <= w_e3_yalmip <= e3_params.d , uncertain(w_e3_yalmip) ]

	    disp('Created YALMIP Constraints.')

	    %Solve Optimization
	    disp('Solving YALMIP Robust Optimization...')
	    ops = sdpsettings('verbose',0);
	    e3_results.sol_yalmip = optimize(l_diag_constr+robust_constr,e3_params.obj_yalmip,ops);

	    if e3_results.sol_yalmip.problem == 0
	    	disp('YALMIP Optimization Solved')
	    else
	    	error('YALMIP Optimization NOT Solved.')
	    end

	    e3_results.Q_yalmip = value(e3_params.Q_yalmip);
	    e3_results.r_yalmip = value(e3_params.r_yalmip);

	    clear Pxw Pxv x_tilde

	end


else
	disp('Experiment not performed.')
end

%% Experiment 5

%Announce Experiment
disp('==========================================================================')
disp('Experiment 5: Comparison of Robust Optimization and Epigraph Reformulation')
disp(' ');


if perform_experiment(5)

	%Load data
	if exist('acc_dsys')
		disp('Data already loaded.')
	else
		load('../data/acc_m1_systems.mat')
	end

	%Defining Constants
	acc_dsys.x0 = [0;1;0];
	acc_dsys.d = 0.1;
	acc_dsys.m = 0.05;
	n = size(acc_dsys.A,1);

	e5_params.n = n;
	e5_params.d = acc_dsys.d;
	e5_params.m = acc_dsys.m;
	e5_params.T = 3;

	% Create Matrices Necessary for Original Skaf and Boyd Method
	%------------------------------------------------------------
	[G_e5,H_e5,Cm_e5,x0m_e5] = create_skaf_n_boyd_matrices(acc_dsys,e5_params.T);
	disp('Created Skaf and Boyd Matrices.')

	% Perform Robust Optimization Using YALMIP "uncertain"
	%-----------------------------------------------------
	e5_params.Q_robust = sdpvar(size(H_e5,2),size(Cm_e5,1),'full');
	e5_params.r_robust = sdpvar(size(H_e5,2),1,'full');

	%Disturbance vectors
	w_e5_yalmip = sdpvar(n*e5_params.T,1,'full');
	v_e5_yalmip = sdpvar(size(acc_dsys.C,1)*e5_params.T,1,'full');
	disp('Created YALMIP Optimization Variables.')

	%Create Expressions Containing Q,r for Optimization
	Pxw = (eye(n*(e5_params.T+1))+H_e5*e5_params.Q_robust*Cm_e5)*G_e5;
	Pxv = H_e5*e5_params.Q_robust;
	x_tilde = (eye(n*(e5_params.T+1)) + H_e5*e5_params.Q_robust*Cm_e5)*x0m_e5 + H_e5*e5_params.r_robust;

	%Create Objective
	R = [zeros(n,n*e5_params.T) eye(n)]; %Create selection matrix
	e5_params.obj_robust = norm( R*(x_tilde + Pxw * w_e5_yalmip + Pxv * v_e5_yalmip) , Inf );
	disp('Created YALMIP Objective.')
	disp(' ')

	%Create Constraints

	%Q is lower diagonal.
	l_diag_constr = [];
	for bl_row_num = 1 : e5_params.T-1
		l_diag_constr = l_diag_constr + [ e5_params.Q_robust(	[(bl_row_num-1)*size(acc_dsys.B,2)+1:bl_row_num*size(acc_dsys.B,2)], ...
																[bl_row_num*size(acc_dsys.C,1)+1:end] ) == 0 ];
    end

    l_diag_constr

    %Robustifying against w and v
    robust_constr = [];
    robust_constr = robust_constr + [ -e5_params.m <= v_e5_yalmip <= e5_params.m , uncertain(v_e5_yalmip) ];
    robust_constr = robust_constr + [ -e5_params.d <= w_e5_yalmip <= e5_params.d , uncertain(w_e5_yalmip) ]

    disp('Created YALMIP Constraints.')

    %Solve Optimization
    disp('Solving YALMIP Robust Optimization...')
    ops = sdpsettings('verbose',0);
    e5_results.sol_robust = optimize(l_diag_constr+robust_constr,e5_params.obj_robust,ops);

    if e5_results.sol_robust.problem == 0
    	disp('YALMIP Robust Optimization Solved')
    else
    	error('YALMIP Robust Optimization NOT Solved.')
    end

    e5_results.Q_robust = value(e5_params.Q_robust);
    e5_results.r_robust = value(e5_params.r_robust);

    clear Pxw Pxv x_tilde

    % Perform Robust Optimization Using the Epigraph Reformulation
    %--------------------------------------------------
    disp(' ')
    disp('Attempting to Perform Worst-Case Optimization using the Epigraph formulation.')
    disp('-----------------------------------------------------------------------------')

    % Create Variables
   	e5_params.Q_epi = sdpvar(size(H_e5,2),size(Cm_e5,1),'full');
	e5_params.r_epi = sdpvar(size(H_e5,2),1,'full');

	%Create Expressions Containing Q,r for Optimization
	Pxw = (eye(n*(e5_params.T+1))+H_e5*e5_params.Q_epi*Cm_e5)*G_e5;
	Pxv = H_e5*e5_params.Q_epi;
	x_tilde = (eye(n*(e5_params.T+1)) + H_e5*e5_params.Q_epi*Cm_e5)*x0m_e5 + H_e5*e5_params.r_epi;


	e5_params.alpha_epi 	= sdpvar(1);
	e5_params.lambda_epi 	= sdpvar(size(R,1),2,'full');
	e5_params.eta_epi 		= sdpvar(size(Pxw,2),2,'full');
	e5_params.gamma_epi 	= sdpvar(size(R,1),2,'full');
	e5_params.theta_epi		= sdpvar(size(Pxv,2),2,'full');
	disp('Created Variables.')

	% Based on the following cvx definitions.
	% variables alpha_e3(1) lambda_e3(size(R,1),2) eta_e3(size(Pxw,2),2) gamma_e3(size(R,1),2) theta_e3(size(Pxv,2),2)

	% Inequality Constraints
	ineq_constrs = [];
	for col_num = 1 : 2
		ineq_constrs = ineq_constrs + [ e5_params.lambda_epi(:,col_num) >= 0 ];
		ineq_constrs = ineq_constrs + [ e5_params.eta_epi(:,col_num) >= 0 ];
		%ineq_constrs = ineq_constrs + [ e5_params.gamma_epi(:,col_num) >= 0 ];
		ineq_constrs = ineq_constrs + [ e5_params.theta_epi(:,col_num) >= 0];
	end

	ineq_constrs = ineq_constrs + ...
					[ -(e5_params.lambda_epi(:,1)-e5_params.lambda_epi(:,2))'*R*x_tilde ...
						- e5_params.d * sum( e5_params.eta_epi(:,1) + e5_params.eta_epi(:,2) ) ...
						- e5_params.m * sum( e5_params.theta_epi(:,1) + e5_params.theta_epi(:,2) ) <= e5_params.alpha_epi ];

	disp('Inequality constraints created.')

	% Equality Constraints.
	eq_constrs = [];
	eq_constrs = eq_constrs + [ sum( e5_params.lambda_epi(:,1) + e5_params.lambda_epi(:,2) ) == 1 ];
	%eq_constrs = eq_constrs + [ e5_params.lambda_epi(:,1)+e5_params.lambda_epi(:,2) == -e5_params.gamma_epi(:,1) + e5_params.gamma_epi(:,2) ];
	eq_constrs = eq_constrs + [ -(e5_params.lambda_epi(:,1)-e5_params.lambda_epi(:,2))'*R * Pxw == (e5_params.eta_epi(:,1) - e5_params.eta_epi(:,2))' ];
	eq_constrs = eq_constrs + [ -(e5_params.lambda_epi(:,1)-e5_params.lambda_epi(:,2))'*R * Pxv == (e5_params.theta_epi(:,1) - e5_params.theta_epi(:,2))' ];
	disp('Equality Constraints Created.')

	%Q is lower diagonal.
	l_diag_constr = [];
	for bl_row_num = 1 : e5_params.T-1
		l_diag_constr = l_diag_constr + [ e5_params.Q_epi(	[(bl_row_num-1)*size(acc_dsys.B,2)+1:bl_row_num*size(acc_dsys.B,2)], ...
																[bl_row_num*size(acc_dsys.C,1)+1:end] ) == 0 ];
    end

    l_diag_constr

	% Optimize
	ops = sdpsettings('verbose',1);
	e5_results.sol_epi = optimize(eq_constrs + ineq_constrs + l_diag_constr,e5_params.alpha_epi,ops)

	% Conclusions
	disp('The epigraph problem does not seem to be solveable in this formulation.')
else
	disp('Experiment not performed.')
end

%% Experiment 6

%Announce Experiment
disp('===============================================================================')
disp('Experiment 6: Comparison of Robust Optimization and Epigraph Reformulation (II)')
disp(' ');


if perform_experiment(6)

	%Load data
	if exist('acc_dsys')
		disp('Data already loaded.')
	else
		load('../data/acc_m1_systems.mat')
	end

	%Defining Constants
	acc_dsys.x0 = [0;1;0];
	acc_dsys.d = 0.1;
	acc_dsys.m = 0.05;
	n = size(acc_dsys.A,1);

	e6_params.n = n;
	e6_params.d = acc_dsys.d;
	e6_params.m = acc_dsys.m;
	e6_params.T = 3;

	% Create Matrices Necessary for Original Skaf and Boyd Method
	%------------------------------------------------------------
	[G_e6,H_e6,Cm_e6,x0m_e6] = create_skaf_n_boyd_matrices(acc_dsys,e6_params.T);
	disp('Created Skaf and Boyd Matrices.')

	% Perform Robust Optimization Using YALMIP "uncertain"
	%-----------------------------------------------------
	e6_params.Q_r = sdpvar(size(H_e6,2),size(Cm_e6,1),'full');
	e6_params.r_r = sdpvar(size(H_e6,2),1,'full');

	%Disturbance vectors
	w_e6 = sdpvar(n*e6_params.T,1,'full');
	v_e6 = sdpvar(size(acc_dsys.C,1)*e6_params.T,1,'full');
	disp('Created YALMIP Optimization Variables.')

	%Create Expressions Containing Q,r for Optimization
	Pxw = (eye(n*(e6_params.T+1))+H_e6*e6_params.Q_r*Cm_e6)*G_e6;
	Pxv = H_e6*e6_params.Q_r;
	x_tilde = (eye(n*(e6_params.T+1)) + H_e6*e6_params.Q_r*Cm_e6)*x0m_e6 + H_e6*e6_params.r_r;

	%Create Objective
	R = [zeros(n,n*e6_params.T) eye(n)]; %Create selection matrix
	e6_params.obj_r = norm( R*(x_tilde + Pxw * w_e6 + Pxv * v_e6) , Inf );
	disp('Created YALMIP Objective.')
	disp(' ')

	%Create Constraints

	%Q is lower diagonal.
	l_diag_constr = [];
	for bl_row_num = 1 : e6_params.T-1
		l_diag_constr = l_diag_constr + [ e6_params.Q_r(	[(bl_row_num-1)*size(acc_dsys.B,2)+1:bl_row_num*size(acc_dsys.B,2)], ...
															[bl_row_num*size(acc_dsys.C,1)+1:end] ) == 0 ];
    end

    % l_diag_constr

    %Robustifying against w and v
    robust_constr = [];
    robust_constr = robust_constr + [ -e6_params.m <= v_e6 <= e6_params.m , uncertain(v_e6) ];
    robust_constr = robust_constr + [ -e6_params.d <= w_e6 <= e6_params.d , uncertain(w_e6) ];

    disp('Created YALMIP Constraints.')

    %Solve Optimization
    disp('Solving YALMIP Robust Optimization...')
    ops = sdpsettings('verbose',0);
    e6_results.sol_robust = optimize(l_diag_constr+robust_constr,e6_params.obj_r,ops);

    if e6_results.sol_robust.problem == 0
    	disp('YALMIP Robust Optimization Solved')
    else
    	error('YALMIP Robust Optimization NOT Solved.')
    end

    e6_results.Q_robust = value(e6_params.Q_r);
    e6_results.r_robust = value(e6_params.r_r);

    disp(' ')

    %clear Pxw Pxv x_tilde w_e6 v_e6

    % Attempt the optimization again with epigraph formulation 
    %---------------------------------------------------------
    disp('Starting Epigraph form of Robust Optimization Problem.')
   	e6_params.alpha = sdpvar(1,1,'full');
    epi_constr = [ e6_params.obj_r <= e6_params.alpha ];

    % Do epigraph minimization
    e6_results.sol_robust_w_epi = optimize(l_diag_constr+robust_constr+epi_constr,e6_params.alpha,ops);

    if e6_results.sol_robust.problem == 0
    	disp('YALMIP Robust Optimization under Epigraph Formulation Solved')
    else
    	error('YALMIP Robust Optimization under Epigraph Formulation NOT Solved.')
    end

    % Save some results.

    e6_results.Q_epi = value(e6_params.Q_r);
    e6_results.r_epi = value(e6_params.r_r);
    e6_results.alpha_epi = value(e6_params.alpha);

    % Compare matrices

    disp(' ')
    disp('Conclusions')

    % disp('e6_results.Q_epi :')
    % e6_results.Q_epi

    % disp('e6_results.Q_robust :')
    % e6_results.Q_robust

    if( sum( sum( e6_results.Q_epi == e6_results.Q_robust ) ) == prod(size(e6_results.Q_epi)) )
    	disp('The Q matrix is the same for both optimizations.')
    else
    	disp('The Q values found in our 2 optimizations are different.')
    end

    if( sum( sum( e6_results.r_epi == e6_results.r_robust ) ) == prod(size(e6_results.r_epi)) )
    	disp('The r matrix is the same for both optimizations.')
    else
    	disp('The r values found in our 2 optimizations are different.')
    end

    disp(['The optimizations yield similar results but we can get an optimal value from alpha (' num2str(value(e6_results.alpha_epi)) ').'])

else
	disp('User decided to skip this experiment.')
end

%% Experiment 7: Varying Horizon

disp('==============================')
disp('Experiment 7: Varying Horizons')
disp(' ')

if perform_experiment(7)

	%Load data
	if exist('acc_dsys')
		disp('Data already loaded.')
	else
		load('../data/acc_m1_systems.mat')
	end

	%Defining Constants
	acc_dsys.x0 = [0;1;0];
	acc_dsys.d = 0.1;
	acc_dsys.m = 0.05;
	n = size(acc_dsys.A,1);

	e7_params.n = n;
	e7_params.d = acc_dsys.d;
	e7_params.m = acc_dsys.m;

	%Change the value of the time horizon T

	e7_params.T_lim = 20;

	for T = 1 : e7_params.T_lim

		% Create Matrices Necessary for Original Skaf and Boyd Method
		%------------------------------------------------------------
		[G_e7,H_e7,Cm_e7,x0m_e7] = create_skaf_n_boyd_matrices(acc_dsys,T);
		disp('Created Skaf and Boyd Matrices.')

		% Perform Robust Optimization Using YALMIP "uncertain"
		%-----------------------------------------------------
		e7_params.Q = sdpvar(size(H_e7,2),size(Cm_e7,1),'full');
		e7_params.r = sdpvar(size(H_e7,2),1,'full');

		%Disturbance vectors
		w_e7 = sdpvar(n*T,1,'full');
		v_e7 = sdpvar(size(acc_dsys.C,1)*T,1,'full');
		disp('Created YALMIP Optimization Variables.')

		%Create Expressions Containing Q,r for Optimization
		Pxw = (eye(n*(T+1))+H_e7*e7_params.Q*Cm_e7)*G_e7;
		Pxv = H_e7*e7_params.Q;
		x_tilde = (eye(n*(T+1)) + H_e7*e7_params.Q*Cm_e7)*x0m_e7 + H_e7*e7_params.r;

		%Create Objective
		R = [zeros(n,n*T) eye(n)]; %Create selection matrix
		e7_params.obj_r = norm( R*(x_tilde + Pxw * w_e7 + Pxv * v_e7) , Inf );
		disp('Created YALMIP Objective.')
		disp(' ')

		%Create Constraints

		%Q is lower diagonal.
		l_diag_constr = [];
		for bl_row_num = 1 : T-1
			l_diag_constr = l_diag_constr + [ e7_params.Q(	[(bl_row_num-1)*size(acc_dsys.B,2)+1:bl_row_num*size(acc_dsys.B,2)], ...
															[bl_row_num*size(acc_dsys.C,1)+1:end] ) == 0 ];
	    end

	    % l_diag_constr

	    %Robustifying against w and v
	    robust_constr = [];
	    robust_constr = robust_constr + [ -e7_params.m <= v_e7 <= e7_params.m , uncertain(v_e7) ];
	    robust_constr = robust_constr + [ -e7_params.d <= w_e7 <= e7_params.d , uncertain(w_e7) ];

	    %Epigraph Constraints
	    e7_params.alpha = sdpvar(1,1,'full');
	    epi_constr = [];
	    epi_constr = [ e7_params.obj_r <= e7_params.alpha ];

	    disp('Created YALMIP Constraints.')

	    %Solve Optimization
	    disp('Solving YALMIP Robust Optimization...')
	    op_num = 0;
	    if T==9
	    	op_num = 2;
	    end
	    ops = sdpsettings('verbose',op_num);
	    e7_results.sol_robust = optimize(l_diag_constr+robust_constr+epi_constr,e7_params.alpha,ops);

	    if e7_results.sol_robust.problem == 0
	    	disp(['YALMIP Robust Optimization #' num2str(T) ' Solved'])
	    else
	    	%error(['YALMIP Robust Optimization #' num2str(T) ' NOT Solved.'])
	    	disp(['YALMIP Robust Optimization #' num2str(T) ' NOT Solved.'])
	    end

	    %Save results
	    e7_results.Q{T} = value(e7_params.Q);
	    e7_results.r{T} = value(e7_params.r);
	    e7_results.alpha{T} = value(e7_params.alpha);

	    clear G_e7 H_e7 Cm_e7 x0m_e7 e7_params.alpha

	end
else
	disp('User decided to skip this experiment.')
end

%% Experiment 8: Comparing Standard Eq. Perf. with the results of Skaf and Boyd

disp('==============================')
disp('Experiment 8: Comparing Design wrt Original Eq. Perf. to New One with Skaf + Boyd')
disp('- Using Error System for Skaf and Boyd Dynamics Now')
disp(' ')

if perform_experiment(8)

	%Load data
	if exist('acc_dsys')
		disp('Data already loaded.')
	else
		load('../data/acc_m1_systems.mat')
	end

	%Defining Constants
	acc_dsys.x0 = [0;1;0];
	acc_dsys.d = 0.1;
	acc_dsys.m = 0.05;
	n = size(acc_dsys.A,1);

	%Create Error System Dynamics
	acc_error_dsys = acc_dsys;

	acc_error_dsys.B = eye(size(acc_dsys.B,1));

	e8_params.n = n;
	e8_params.d = acc_dsys.d;
	e8_params.m = acc_dsys.m;

	e8_params.perf_level = 1;

	%We are only concerned with the next step.
	e8_params.T = 1;

	% Create Matrices Necessary for Original Skaf and Boyd Method
	%------------------------------------------------------------
	[G_e8,H_e8,Cm_e8,x0m_e8] = create_skaf_n_boyd_matrices(acc_error_dsys,e8_params.T);
	disp('Created Skaf and Boyd Matrices.')

	% Perform Robust Optimization Using YALMIP "uncertain"
	%-----------------------------------------------------
	e8_params.Q = sdpvar(size(H_e8,2),size(Cm_e8,1),'full');
	e8_params.r = sdpvar(size(H_e8,2),1,'full');

	%Disturbance vectors
	w_e8 = sdpvar(n*e8_params.T,1,'full');
	v_e8 = sdpvar(size(acc_error_dsys.C,1)*e8_params.T,1,'full');
	disp('Created YALMIP Optimization Variables.')

	%Create Expressions Containing Q,r for Optimization
	Pxw = (eye(n*(e8_params.T+1))+H_e8*e8_params.Q*Cm_e8)*G_e8;
	Pxv = H_e8*e8_params.Q;
	x_tilde = (eye(n*(e8_params.T+1)) + H_e8*e8_params.Q*Cm_e8)*x0m_e8 + H_e8*e8_params.r;

	%Create Objective
	R = [zeros(n,n*e8_params.T) eye(n)]; %Create selection matrix
	e8_params.obj_r = norm( R*(x_tilde + Pxw * w_e8 + Pxv * v_e8) , Inf );
	disp('Created YALMIP Objective.')
	disp(' ')

	%Create Constraints

	%Q is lower diagonal.
	l_diag_constr = [];
	for bl_row_num = 1 : e8_params.T-1
		l_diag_constr = l_diag_constr + [ e8_params.Q(	[(bl_row_num-1)*size(acc_error_dsys.B,2)+1:bl_row_num*size(acc_error_dsys.B,2)], ...
														[bl_row_num*size(acc_error_dsys.C,1)+1:end] ) == 0 ];
    end

    % l_diag_constr

    %Robustifying against w and v
    robust_constr = [];
    robust_constr = robust_constr + [ -e8_params.m <= v_e8 <= e8_params.m , uncertain(v_e8) ];
    robust_constr = robust_constr + [ -e8_params.d <= w_e8 <= e8_params.d , uncertain(w_e8) ];

    %Epigraph Constraints
    e8_params.alpha = sdpvar(1,1,'full');
    epi_constr = [];
    epi_constr = [ e8_params.obj_r <= e8_params.alpha ];

    disp('Created YALMIP Constraints.')

    %Solve Optimization
    disp('Solving YALMIP Robust Optimization...')
    op_num = 0;
    ops = sdpsettings('verbose',op_num);
    e8_results.sol_robust = optimize(l_diag_constr+robust_constr+epi_constr,e8_params.alpha,ops);

    if e8_results.sol_robust.problem == 0
    	disp(['YALMIP Robust Optimization Solved'])
    else
    	%error(['YALMIP Robust Optimization #' num2str(T) ' NOT Solved.'])
    	disp(['YALMIP Robust Optimization NOT Solved.'])
    end

    %Save results
    e8_results.Q = value(e8_params.Q);
    e8_results.r = value(e8_params.r);
    e8_results.alpha = value(e8_params.alpha);

    clear G_e8 H_e8 Cm_e8 x0m_e8 e8_params.alpha

    % Can the above be compacted into a single function?
    %---------------------------------------------------

    disp(' ')
    disp('Attempting to achieve the same result using a function...')

    e8b_results = generate_skaf_controller(acc_error_dsys,e8_params.T,0);

    disp('Function works!')
    disp('+++++++++++++++')
    disp(' ')

    % Perform Design with the Previous Method
    %----------------------------------------

    e8_params.L = sdpvar(size(acc_error_dsys.A,1),size(acc_error_dsys.C,1),'full');
    e8_params.old_obj = norm( acc_error_dsys.A - e8_params.L*acc_error_dsys.C , Inf ) + (e8_params.m/e8_params.perf_level)*norm(e8_params.L,Inf);
    e8_results.std_eq_perf.sol = optimize([],e8_params.old_obj,ops);

    e8_results.std_eq_perf.L = value(e8_params.L);
    e8_results.std_eq_perf.opt_obj = value(e8_params.old_obj);

    % Evaluate the two feedback laws
    %-------------------------------

    disp('Compare our 2 feedback laws.')
    disp('Which one truly allows for the best operation given worst case noise?')
    disp(' ')

    fb_gains(:,:,1) = e8b_results.F;
    fb_gains(:,:,2) = e8_results.std_eq_perf.L;

    fb_affine_terms(:,:,1) = e8b_results.u0;
    fb_affine_terms(:,:,2) = 0;

    %Create sdpvars for noise
    e8_params.w = sdpvar(size(acc_error_dsys.A,1),1,'full');
    e8_params.v = sdpvar(size(acc_error_dsys.C,1),1,'full');

    for fb_num = 1 : size(fb_gains,3)
 		%Create Objective
 		comparison_obj = norm( acc_error_dsys.A*acc_error_dsys.x0 + e8_params.w + ...
 								fb_gains(:,:,fb_num) * acc_dsys.C * acc_error_dsys.x0 + ...
 								fb_gains(:,:,fb_num)*e8_params.v + fb_affine_terms(:,:,fb_num) , Inf );

 		%Create robust constraints

 		alpha1 = sdpvar(1,1,'full');

	    %Robustifying against w and v
	    robust_constr = [];
	    robust_constr = robust_constr + [ -e8_params.m <= e8_params.v <= e8_params.m , uncertain(e8_params.v) ];
	    robust_constr = robust_constr + [ -e8_params.d <= e8_params.w <= e8_params.d , uncertain(e8_params.w) ];

	    e8_results.comp_fbs.sol(fb_num) = optimize(robust_constr + [ comparison_obj <= alpha1 ],alpha1,sdpsettings('verbose',2));
	    e8_results.comp_fbs.opt_obj(fb_num) = value(alpha1);

 	end

 	disp('Showing the worst-case values of the error at the next step for each method.')
 	e8_results.comp_fbs.opt_obj

else
	disp('User decided to skip this experiment.')
end

%% Experiment 9: Varying Horizon w/ proper model structure

disp('==========================================================')
disp('Experiment 9: Varying Horizon for Observer Feedback Design')
disp('- Using Error System for Skaf and Boyd Dynamics')
disp('- Function used')

if perform_experiment(9)

	%Defining Constants
	acc_dsys.x0 = [0;1;0];
	acc_dsys.d = 0.1;
	acc_dsys.m = 0.05;
	n = size(acc_dsys.A,1);

	%Create Error System Dynamics
	acc_error_dsys = acc_dsys;

	acc_error_dsys.B = eye(size(acc_dsys.B,1));

	e9_params.n = n;
	e9_params.d = acc_dsys.d;
	e9_params.m = acc_dsys.m;
	e9_params.perf_level = 1;

	e9_params.T = 20;

	for T = 1 : e9_params.T
		%Design Controller
		e9_results.skaf_results(T) = generate_skaf_controller(acc_error_dsys,T,0);
	
		e9_results.obj_vals(T) = e9_results.skaf_results(T).opt_obj;
		if mod(T,e9_params.T/20) == 0
			disp([ num2str(T) ' iterations passed!'])
		end
	end

	%Plot Results
	figure;
	stem(e9_results.obj_vals)

	save('data/e9_results.mat','e9_results')

else
	disp('User decided to skip this experiment.')
end

%% Experiment 10: Making x0 a Robust Variable in optimization

disp('=====================================')
disp('Experiment 10: Robustifying w.r.t. x0')
disp('- Using Error System for Skaf and Boyd Dynamics')
disp('- Making x0 also an optimization variable')
disp(' ')

if perform_experiment(10)

	%Defining Constants
	acc_dsys.x0 = sdpvar(3,1,'full');
	acc_dsys.d = 0.1;
	acc_dsys.m = 0.05;
	n = size(acc_dsys.A,1);

	%Create Error System Dynamics
	acc_error_dsys = acc_dsys;

	acc_error_dsys.B = eye(size(acc_dsys.B,1));

	e10_params.n = n;
	e10_params.d = acc_dsys.d;
	e10_params.m = acc_dsys.m;
	e10_params.perf_level = 1;

	e10_params.T = 1;

	% Create Matrices Necessary for Original Skaf and Boyd Method
	%------------------------------------------------------------
	[G_e10,H_e10,Cm_e10,x0m_e10] = create_skaf_n_boyd_matrices(acc_error_dsys,e10_params.T);
	disp('Created Skaf and Boyd Matrices.')

	% Perform Robust Optimization Using YALMIP "uncertain"
	%-----------------------------------------------------
	e10_params.Q = sdpvar(size(H_e10,2),size(Cm_e10,1),'full');
	e10_params.r = sdpvar(size(H_e10,2),1,'full');

	%Disturbance vectors
	w_e10 = sdpvar(n*e10_params.T,1,'full');
	v_e10 = sdpvar(size(acc_error_dsys.C,1)*e10_params.T,1,'full');
	disp('Created YALMIP Optimization Variables.')

	%Create Expressions Containing Q,r for Optimization
	Pxw = (eye(n*(e10_params.T+1))+H_e10*e10_params.Q*Cm_e10)*G_e10;
	Pxv = H_e10*e10_params.Q;
	x_tilde = (eye(n*(e10_params.T+1)) + H_e10*e10_params.Q*Cm_e10)*x0m_e10 + H_e10*e10_params.r;

	%Create Objective
	R = [zeros(n,n*e10_params.T) eye(n)]; %Create selection matrix
	e10_params.obj_r = norm( R*(x_tilde + Pxw * w_e10 + Pxv * v_e10) , Inf );
	disp('Created YALMIP Objective.')
	disp(' ')

	%Create Constraints

	%Q is lower diagonal.
	l_diag_constr = [];
	for bl_row_num = 1 : e10_params.T-1
		l_diag_constr = l_diag_constr + [ e10_params.Q(	[(bl_row_num-1)*size(acc_error_dsys.B,2)+1:bl_row_num*size(acc_error_dsys.B,2)], ...
														[bl_row_num*size(acc_error_dsys.C,1)+1:end] ) == 0 ];
    end

    % l_diag_constr

    %Robustifying against w and v
    robust_constr = [];
    robust_constr = robust_constr + [ -e10_params.m <= v_e10 <= e10_params.m , uncertain(v_e10) ];
    robust_constr = robust_constr + [ -e10_params.d <= w_e10 <= e10_params.d , uncertain(w_e10) ];
    robust_constr = robust_constr + [ -e10_params.perf_level <= acc_error_dsys.x0 <= e10_params.perf_level , uncertain(acc_error_dsys.x0)];

    %Epigraph Constraints
    e10_params.alpha = sdpvar(1,1,'full');
    epi_constr = [];
    epi_constr = [ e10_params.obj_r <= e10_params.alpha ];

    disp('Created YALMIP Constraints.')

    %Solve Optimization
    disp('Solving YALMIP Robust Optimization...')
    op_num = 0;
    ops = sdpsettings('verbose',op_num);
    e10_results.sol_robust = optimize(l_diag_constr+robust_constr+epi_constr,e10_params.alpha,ops);

    if e10_results.sol_robust.problem == 0
    	disp(['YALMIP Robust Optimization Solved'])
    else
    	error(['YALMIP Robust Optimization #' num2str(T) ' NOT Solved.'])
    	%disp(['YALMIP Robust Optimization NOT Solved.'])
    end

    %Save results
    e10_results.Q = value(e10_params.Q);
    e10_results.r = value(e10_params.r);
    e10_results.alpha = value(e10_params.alpha);

    e10_results.F = value( (pinv(value(eye(size(e10_params.Q,1)) + e10_params.Q*Cm_e10*H_e10)) ) * e10_params.Q);
	e10_results.u0 = value((eye(size(e10_results.F,1)) + e10_results.F*Cm_e10*H_e10) * e10_params.r);

	save('data/afhc_e10.mat','e10_results')

	%Display something about the results.
	disp('Showing the values of ||.||_infty for this method.')
	disp([ num2str(e10_results.alpha)])
else
	disp('User decided to skip this experiment.')
end

%% Experiment 10: Making x0 a Robust Variable in optimization

disp('================================================')
disp('Experiment 11: Robustifying w.r.t. x0, Varying T')
disp('- Using Error System for Skaf and Boyd Dynamics')
disp('- Making x0 also an optimization variable (for robust optim.)')
disp(' ')

experim_num = 11;

if perform_experiment(experim_num)

	%Defining Constants
	acc_dsys.x0 = sdpvar(3,1,'full');
	acc_dsys.d = 0.1;
	acc_dsys.m = 0.05;
	n = size(acc_dsys.A,1);

	%Create Error System Dynamics
	acc_error_dsys = acc_dsys;

	acc_error_dsys.B = eye(size(acc_dsys.B,1));

	exp_params{experim_num}.n = n;
	exp_params{experim_num}.d = acc_dsys.d;
	exp_params{experim_num}.m = acc_dsys.m;
	exp_params{experim_num}.perf_level = 1;

	exp_params{experim_num}.T_lim = 8;

	for T = 1: exp_params{experim_num}.T_lim

		exp_params{experim_num}.T = T;

		% Create Matrices Necessary for Original Skaf and Boyd Method
		%------------------------------------------------------------
		[G{experim_num},H{experim_num},Cm{experim_num},x0m{experim_num}] = create_skaf_n_boyd_matrices(acc_error_dsys,exp_params{experim_num}.T);
		disp('Created Skaf and Boyd Matrices.')

		% Perform Robust Optimization Using YALMIP "uncertain"
		%-----------------------------------------------------
		exp_params{experim_num}.Q = sdpvar(size(H{experim_num},2),size(Cm{experim_num},1),'full');
		exp_params{experim_num}.r = sdpvar(size(H{experim_num},2),1,'full');

		%Disturbance vectors
		w{experim_num} = sdpvar(n*exp_params{experim_num}.T,1,'full');
		v{experim_num} = sdpvar(size(acc_error_dsys.C,1)*exp_params{experim_num}.T,1,'full');
		disp('Created YALMIP Optimization Variables.')

		%Create Expressions Containing Q,r for Optimization
		Pxw = (eye(n*(exp_params{experim_num}.T+1))+H{experim_num}*exp_params{experim_num}.Q*Cm{experim_num})*G{experim_num};
		Pxv = H{experim_num}*exp_params{experim_num}.Q;
		x_tilde = (eye(n*(exp_params{experim_num}.T+1)) + H{experim_num}*exp_params{experim_num}.Q*Cm{experim_num})*x0m{experim_num} + H{experim_num}*exp_params{experim_num}.r;

		%Create Objective
		R = [zeros(n,n*exp_params{experim_num}.T) eye(n)]; %Create selection matrix
		exp_params{experim_num}.obj_r = norm( R*(x_tilde + Pxw * w{experim_num} + Pxv * v{experim_num}) , Inf );
		disp('Created YALMIP Objective.')
		disp(' ')

		%Create Constraints

		%Q is lower diagonal.
		l_diag_constr = [];
		for bl_row_num = 1 : exp_params{experim_num}.T-1
			l_diag_constr = l_diag_constr + [ exp_params{experim_num}.Q(	[(bl_row_num-1)*size(acc_error_dsys.B,2)+1:bl_row_num*size(acc_error_dsys.B,2)], ...
																			[bl_row_num*size(acc_error_dsys.C,1)+1:end] ) == 0 ];
	    end

	    % l_diag_constr

	    %Robustifying against w and v
	    robust_constr = [];
	    robust_constr = robust_constr + [ -exp_params{experim_num}.m 			<= v{experim_num} 	<= exp_params{experim_num}.m , uncertain(v{experim_num}) ];
	    robust_constr = robust_constr + [ -exp_params{experim_num}.d 			<= w{experim_num} 	<= exp_params{experim_num}.d , uncertain(w{experim_num}) ];
	    robust_constr = robust_constr + [ -exp_params{experim_num}.perf_level 	<= acc_error_dsys.x0 <= exp_params{experim_num}.perf_level , uncertain(acc_error_dsys.x0)];

	    %Epigraph Constraints
	    exp_params{experim_num}.alpha = sdpvar(1,1,'full');
	    epi_constr = [];
	    epi_constr = [ exp_params{experim_num}.obj_r <= exp_params{experim_num}.alpha ];

	    disp('Created YALMIP Constraints.')

	    %Solve Optimization
	    disp('Solving YALMIP Robust Optimization...')
	    op_num = 2;
	    ops = sdpsettings('verbose',op_num);
	    exp_results{experim_num}.sol_robust{T} = optimize(l_diag_constr+robust_constr+epi_constr,exp_params{experim_num}.alpha,ops);

	    if exp_results{experim_num}.sol_robust{T}.problem == 0
	    	disp(['YALMIP Robust Optimization Solved'])
	    else
	    	%error(['YALMIP Robust Optimization #' num2str(T) ' NOT Solved.'])
	    	error(['YALMIP Robust Optimization NOT Solved (' yalmiperror(exp_results{experim_num}.sol_robust{T}.problem) ').' ])
	    end

	    %Save results
	    exp_results{experim_num}.Q{T} = value(exp_params{experim_num}.Q);
	    exp_results{experim_num}.r{T} = value(exp_params{experim_num}.r);
	    exp_results{experim_num}.opt_obj(T) = value(exp_params{experim_num}.alpha);

	    exp_results{experim_num}.F{T} = value( (pinv(value(eye(size(exp_params{experim_num}.Q,1)) + exp_params{experim_num}.Q*Cm{experim_num}*H{experim_num})) ) * exp_params{experim_num}.Q);
		exp_results{experim_num}.u0{T} = value((eye(size(exp_results{experim_num}.F{T},1)) + exp_results{experim_num}.F{T}*Cm{experim_num}*H{experim_num}) * exp_params{experim_num}.r);


		% Compare to function output
		exp_results{experim_num}.fcn{T} = generate_skaf_controller( acc_error_dsys , T , 0 , exp_params{experim_num}.perf_level );
		exp_results{experim_num}.fcn_opt_objs(T) = exp_results{experim_num}.fcn{T}.opt_obj;
	end

	save(['data/afhc_e' num2str(experim_num) '.mat'],'exp_results')

	%Display something about the results.
	% disp('Showing the values of ||.||_infty for this method.')
	% disp([ num2str(exp_results{11}.alpha)])

	figure;
	plot(exp_results{experim_num}.opt_obj)
	xlabel('Time Horizon T')
	ylabel('Optimal Objective Value (from Skaf Design)')
	title('Skaf Design''s Objective Value as Horizon Changes')
	axis([1 8 0 exp_results{experim_num}.opt_obj(8)+1])

	figure;
	plot(exp_results{experim_num}.opt_obj)
	xlabel('Time Horizon T')
	ylabel('Optimal Objective Value (from Skaf Design)')
	title('Skaf Design''s Objective Value as Horizon Changes')
	axis([1 8 0 exp_results{experim_num}.opt_obj(8)+1])

else
	disp('User decided to skip this experiment.')
end

%% Remove Functions
rmpath('./functions')

