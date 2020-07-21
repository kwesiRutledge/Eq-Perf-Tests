function [results] = observer_comparison81( varargin )
	%observer_comparison81.m
	%Description:
	%	Playing around with the 
	%

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	test_name = 'observer_comparison81';

	disp(['Beginning ' test_name '.' ])
	disp(' ')

	disp('Defining Constants.')

	%A = [1,0.5;0,0.8];
	A = eye(2);
	n = size(A,1);
	
	B = eye(2);
	m = size(B,2);

	f = zeros(n,1);
	C = eye(n);

	eta_w = 0.4;
	Pw_template = Polyhedron('lb',-eta_w*ones(1,n),'ub',eta_w*ones(1,n));
	Pw1 = Pw_template + [eta_w;0];
	Pw2 = Pw_template + [0;eta_w];
	Pv = Polyhedron('lb',zeros(1,n),'ub',zeros(1,n));

	ad1 = Aff_Dyn(A,B,f,C,Pw1,Pv);
	ad2 = Aff_Dyn(A,B,f,C,Pw2,Pv);

	sys = LCSAS( [ad1,ad2] , Language([1,1,1],[2,2,2]) );

	eta_u = 1;
	Pu = Polyhedron('lb',-eta_u*ones(1,m),'ub',eta_u*ones(1,m));

	eta_x0 = 0.1;
	Px0 = Polyhedron('lb',-eta_x0*ones(1,n),'ub',eta_x0*ones(1,n));

	verbosity = 1;

	flags.ConstructBeliefGraph = true;

    % Save Results
    results.Parameters.System = sys;
    results.Parameters.Verbosity = verbosity;
    results.Parameters.flags = flags;
    results.Parameters.Pu = Pu;
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Display Belief Graph %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%

	if flags.ConstructBeliefGraph
		%% Get the Belief Graph for this Simple 2D System
		BG = BeliefGraph(sys,Pu,Px0,'fb_method','state','verbosity',verbosity);

		figure;
		BG.plot()
		results.Experiment1.BG = BG;
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Optimize Using Binary Variables to Incorporate Reachability/Consistency Results %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	disp('Beginning test 1.')
	disp('Creating an optimization that uses the mixed integer linear variables.')

	%% Constants %%

	plane_pos = 2.0;
	P_target1 = Polyhedron('A',[-1,0],'b',-plane_pos);
	P_target2 = Polyhedron('A',[0,-1],'b',-plane_pos);

	n_u = size(sys.Dyn(1).B,2);
	n_x = size(sys.Dyn(1).A,1);
	n_y = size(sys.Dyn(1).C,1);
	n_w = size(sys.Dyn(1).B_w,2);
	n_v = size(sys.Dyn(1).C_v,2);

	constraint_generator = constr_gen();

	ops = sdpsettings('verbose',verbosity);

	sigma0 = sys.L.words{1};
	T = length(sigma0);

	[H,S,C_bar,J,f_bar] = sys.get_mpc_matrices('word',sigma0);

	P_eta = 1;
	P_wT = 1;
	for symb_idx = 1:length(sigma0)
		%Append the appropriate P_w values
		P_eta = P_eta * sys.Dyn( sigma0(symb_idx) ).P_w;
		P_wT = P_wT * sys.Dyn( sigma0(symb_idx) ).P_w;
	end
	P_eta = P_eta * Px0;

	PuT = 1;
	for symb_idx = 1:length(sigma0)
		%Append the appropriate P_w values
		PuT = PuT*Pu;
	end

	select_m = @(t,T_r) [zeros(n_x,t*n_x), eye(n_x), zeros(n_x,(T_r-t)*n_x) ];

	M = (10^(3)) * ones(size(P_target1.b));

	%% Optimization Variables %%

	Q = sdpvar(n_u*T,n_x*(T+1),'full');
	r = sdpvar(n_u*T,1,'full');

	target1_flag = binvar(1);
	target2_flag = binvar(1);

	Pi1 = sdpvar(size(P_target1.A,1),size(P_eta.A,1),'full');
	Pi2 = sdpvar(size(P_target2.A,1),size(P_eta.A,1),'full');

	Piu1 = sdpvar(size(PuT.A,1),size(P_eta.A,1),'full');
	Piu2 = sdpvar(size(PuT.A,1),size(P_eta.A,1),'full');

	%% Create Constraints %%

	binvar_constr = [ target1_flag + target2_flag == 1 ];

	positivity_constr = [ Pi1 >= 0 , Pi2 >= 0 , Piu1 >= 0 , Piu2 >= 0 ];

	tempI = eye(size(Q,2));
	containment_constr = ...
		[ Pi1*P_eta.A == P_target1.A*select_m(T,T)*[ (tempI+S*Q)*H , (tempI+S*Q)*J ] ] + ...
		[ Pi1*P_eta.b <= P_target1.b + (1-target1_flag)*M - P_target1.A*select_m(T,T)*S*r ] + ...
		[ Pi2*P_eta.A == P_target2.A*select_m(T,T)*[ (tempI+S*Q)*H , (tempI+S*Q)*J ] ] + ...
		[ Pi2*P_eta.b <= P_target2.b + (1-target2_flag) * M - P_target2.A*select_m(T,T)*S*r ];
		
	gain_constrs = [];
	for tau = 0:T-1
		gain_constrs = gain_constrs + [ Q(n_u*tau+[1:n_u],[n_x*(tau+1)+1:end]) == 0 ];
	end

	% control_input_constr = ...
	% 	[ Piu1*P_eta.A == PuT.A*[ Q*H , Q*J ] ] + ...
	% 	[ Piu1*P_eta.b <= PuT.b + target1_flag * M - PuT.A*r ] + ...
	% 	[ Piu2*P_eta.A == PuT.A*[ Q*H , Q*J ] ] + ...
	% 	[ Piu2*P_eta.b <= PuT.b + target2_flag * M - PuT.A*r ];

	control_input_constr = ...
		[ Piu1*P_eta.A == PuT.A*[ Q*H , Q*J ] ] + ...
		[ Piu1*P_eta.b <= PuT.b - PuT.A*r ] + ...
		[ Piu2*P_eta.A == PuT.A*[ Q*H , Q*J ] ] + ...
		[ Piu2*P_eta.b <= PuT.b - PuT.A*r ];

	%% Synthesize Gains %%

	%target1_flag should be assigned 1 and target2_flag should be assigned 0

	opt_out = optimize(	binvar_constr + positivity_constr + containment_constr + gain_constrs + control_input_constr, ...
						[] , ...
						ops);

	%Extract Gains
	Q_opt = value(Q);
	r_opt = value(r);

	tempI = eye(size(Q_opt,1));
	K = (tempI+Q_opt*S)^(-1) * Q_opt;
	k = (tempI+Q_opt*S)^(-1)*r_opt;

	%Kp = [ K, zeros(size(K,1),n_x) ];

	tempI = eye(size(S,1));
	x = (H + S*K*(tempI-S*K)^(-1)*H)*P_wT + ( eye(size(J,1)) + S*K*(tempI-S*K)^(-1) ) * (J*Px0 + S*k);
	x.computeVRep();

	results.Experiment2.opt_out = opt_out;
	results.Experiment2.K = K;
	results.Experiment2.k = k;


	disp(['value(target1_flag) = ' num2str(value(target1_flag)) ])
	disp(['value(target2_flag) = ' num2str(value(target2_flag)) ])



	%% Plotting Results %%

	eta_domain = 10;
	Domain = Polyhedron('lb',-eta_domain*ones(1,n_x),'ub',eta_domain*ones(1,n_x));

	figure;
	hold on;
	plot(Domain.intersect(P_target1),'Alpha',0.5)
	plot(Domain.intersect(P_target2),'Alpha',0.5)

	plot(select_m(T,T)*x,'color','cyan')

	% axis([-4 4 -4 4])

    %% Compute Gain Synthesis for this Problem when Target Region is a Union %%
    
    %Clear some old variables out
    clear H S C_bar J f_bar P_eta P_wT Q r Pi1 Pi2 Piu1 Piu2 target_flag positivity_constr K k
    
    P_targets = [ P_target1 , P_target2 ];
	
    %T is already defined
    
    %Create Constants
    for word_idx = 1:length(sys.L.words)
       temp_word = sys.L.words{word_idx};
       
       % Get MPC Matrices
       [H{word_idx},S{word_idx},C_bar{word_idx},J{word_idx},f_bar{word_idx}] = sys.get_mpc_matrices('word',temp_word);
       
       %Create disturbance profile using temp_word
       temp_P_wT = 1;
       temp_P_eta = 1;
       for symb_idx = 1:length(temp_word)
           %Append the appropriate P_w values
           temp_P_eta = temp_P_eta * sys.Dyn( temp_word(symb_idx) ).P_w;
           temp_P_wT = temp_P_wT * sys.Dyn( temp_word(symb_idx) ).P_w;
       end
       temp_P_eta = temp_P_eta * Px0;
       
       P_eta{word_idx} = temp_P_eta;
       P_wT{word_idx} = temp_P_wT;
       
    end
    
    %P_uT already exists
    
    M = (10^(3)) * ones(size(P_targets(1).b));
    
    %% Optimization Variables %%

    %Gain variables are needed for each word
    for word_idx = 1:sys.L.cardinality()
        Q{word_idx} = sdpvar(n_u*T,n_x*(T+1),'full');
        r{word_idx} = sdpvar(n_u*T,1,'full');
        
        for target_idx = 1:length(P_targets)
            %for each target associate a dual variable for inclusion
            Pi1{word_idx,target_idx} = sdpvar(size(P_targets(target_idx).A,1),size(P_eta{word_idx}.A,1),'full');
%             Pi2{word_idx,target_idx} = sdpvar(size(P_targets(2).A,1),size(P_eta.A,1),'full');
        end
        
        %For each target there should be a boolean variable
        target_flag{word_idx} = binvar(length(P_targets),1);
        
        %Create dual variables for input constraints
        Piu{word_idx} = sdpvar(size(PuT.A,1),size(P_eta{word_idx}.A,1),'full');
    end
    
    %% Create Constraints %%
    
    %Create placeholder constraints
    constr_binvar_reach_atarget = [];
    positivity_constr = [];
    containment_constr = [];
    gain_constrs = [];
    control_input_constr = [];
    
    %For each gain, make sure:
    for word_idx = 1:sys.L.cardinality()
        %Make sure that one of the targets is reached for each gain.
        constr_binvar_reach_atarget = constr_binvar_reach_atarget + [ sum( target_flag{word_idx} ) == 1 ];
        
        %Make sure that the dual variables are positive.
        for target_idx = 1:length(P_targets)
            positivity_constr = positivity_constr + [ Pi1{word_idx,target_idx} >= 0 ]; 
        end
        positivity_constr = positivity_constr + [ Piu{word_idx} >= 0 ];
        
        %Write down the containment constraints in terms of the dual
        %variables.
        tempI = eye(size(Q{word_idx},2));
        Sw = S{word_idx};
        Hw = H{word_idx};
        Jw = J{word_idx};
        for target_idx = 1:length(P_targets)
            containment_constr = containment_constr + ...
                [ Pi1{word_idx,target_idx}*P_eta{word_idx}.A == P_targets(target_idx).A*select_m(T,T)*[ (tempI+Sw*Q{word_idx})*Hw , (tempI+Sw*Q{word_idx})*Jw ] ] + ...
                [ Pi1{word_idx,target_idx}*P_eta{word_idx}.b <= P_targets(target_idx).b + (1-target_flag{word_idx}(target_idx))*M - P_targets(target_idx).A*select_m(T,T)*Sw*r{word_idx} ];
        end
        
        %Constrain the causality of the controller
        for tau = 0:T-1
            gain_constrs = gain_constrs + [ Q{word_idx}(n_u*tau+[1:n_u],[n_x*(tau+1)+1:end]) == 0 ];
        end

        %Constrain the input size to always be within the given bounds
        %(PuT)
        control_input_constr = control_input_constr + ...
            [ Piu{word_idx}*P_eta{word_idx}.A == PuT.A*[ Q{word_idx}*Hw , Q{word_idx}*Jw ] ] + ...
            [ Piu{word_idx}*P_eta{word_idx}.b <= PuT.b - PuT.A*r{word_idx} ];
        
    end
    
    %% Synthesize Gains for Experiment 2 %%
    opt_out = optimize(	constr_binvar_reach_atarget + positivity_constr + containment_constr + gain_constrs + control_input_constr, ...
						[] , ...
						ops);
    
    results.Experiment2.opt_out = opt_out;
                    
    %Extract Gains
    for word_idx = 1:sys.L.cardinality()
        Qw = value(Q{word_idx});
        rw = value(r{word_idx});

        Sw = S{word_idx};
        
        tempI = eye(size(Qw,1));
        Kw = (tempI+Qw*Sw)^(-1) * Qw;
        kw = (tempI+Qw*Sw)^(-1) * rw;

        K{word_idx} = Kw;
        k{word_idx} = kw;
        
    end
    
    results.Experiment2.K = K;
    results.Experiment2.k = k;
    
    %Extract Binary choices
    for word_idx = 1:sys.L.cardinality()
        target_choice_w = value(target_flag{word_idx});
        target_choice{word_idx} = target_choice_w;
    end
    results.Experiment2.TargetChoices = target_choice;
        
end