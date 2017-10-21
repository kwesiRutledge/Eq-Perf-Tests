% observer_comparison3.m
%
%	Trying to show that 

%% Constants

%Creating a system matrix with arbitrary roots.

a = [ 1 3 3 ]

A = [ 0    1    0 ;
	  0    0    1 ;
	  a(3) a(2) a(1) ];

C{1} = [1 0 0];
C{2} = [0 1 0];
C{3} = [0 0 1];

E = [1 0 0]';

%Save this example system into a struct
ex1.A = A;
ex1.B = eye(size(A,1)); 	%We don't need a B matrix here, but I am applying it so that the function doesn't throw any errors.
ex1.E = E;
ex1.m = 0;
ex1.d = 0.1;
ex1.x0 = Inf;

pl = 1;

%% Experiment 1

disp('==============')
disp('Experiment 3.1')
disp('What is the effect of the observed state C on the optimal value of ''guaruntee''-able error?')

%Looping through different C matrices
for c_num = 1:3

	%sol_skaf = 

	%Looping through different time horizons.
	for t =  1 : 10

		%Define C
		ex1.C = C{c_num};

		%Calculate Skaf and Yong Solutions at this time step
		sol_skaf = generate_skaf_controller(ex1,t,verbosity,'PL',pl);
		%sol_yong = generate_yong_controller(ex1,t,verbosity,'PL',pl);

		%Save data to vector
		skaf_opt(t) = sol_skaf.opt_obj;
		%yong_opt(t) = sol_yong.opt_obj;

		%Triangle Inequality
		L 		= sdpvar(size(ex1.A,1),size(ex1.C,1),'full');
		alpha_0 = sdpvar(1,1,'full');

		%Create optimization objective.
		e_t = norm( (ex1.A-L*ex1.C) , Inf )^t;
		for k = 0:t-1
			e_t = e_t + norm( ( (ex1.A-L*ex1.C)^k ) * L , Inf ) * (ex1.m/pl) + norm( ( (ex1.A-L*ex1.C)^k ) * acc.E , Inf )*(ex1.d/pl) ;
		end
		ti_sol2(t) = optimize([],e_t,sdpsettings('verbose',verbosity));

		tri_ineq_solved(t) = (ti_sol2(t).problem==0);
		tri_ineq_opt(t) = value(e_t);

		disp( [ num2str(t) ' iterations finished.' ])

	end

	exp_results{c_num} = skaf_opt;

	% Plot results

	figure;
	hold on;
	plot(skaf_opt)
	%plot(yong_opt)
	%plot(tri_ineq_opt)
	plot(ones(1,10),':')

	title('Guarunteed Error at Time Step T')
	xlabel('Time Horizon (T)')
	ylabel('Error Magnitude $||e(t)||_{\infty}$','Interpreter','latex')
	legend('Skaf','Desired Error')
	%legend('Skaf','Yong','Triangle Inequality','Desired Error')
		
end

%Plot all 3 together
symb_list = {'o','b','x'};
figure;
hold on;
for c_num = 1:length(exp_results)
	plot(exp_results{c_num},symb_list{c_num})
end
plot(ones(1,10),':')

title('Effect of Different Measurement Matrices on Guarunteed Error')
xlabel('Time Horizon (T)')
ylabel('Error Magnitude $||e(t)||_{\infty}$','Interpreter','latex')
legend('C1','C2','C3','Desired Error')
