function [results] = observer_comparison33(varargin)
	%observer_comparison33.m
	%	This script is meant to observe the minimum value of M_3 as M_1 varies
	%	for various problems. For some problems, equalized performance is
	%	possible and for others it is not.

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	load('data/system_examples/acc_p.mat');

	%Create Aff_Dyn object with the data from acc_e
	acc_e = acc;
	acc_e.B = eye(size(acc.A,1));

	acc_ad = Aff_Dyn(acc_e.A,acc_e.B,zeros(size(acc_e.A,1),1), acc_e.C, acc_e.d , acc_e.m, acc_e.E , eye(size(acc.C,1)) );

	L1 = [1,0,1,0,0];

	attempt_these_tests = [1 2 3 4];

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Iterate through Some Values of M1 %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	M1_vals = [0.0:0.01:0.5]; 
	M3_vals = [];

	if any(attempt_these_tests == 1)

		for M1 = M1_vals

			[opt_out,estima] = acc_ad.free_rec_design_pb('Min_M3',M1,M1+2,L1,'verbosity',0);

			if opt_out.problem == 0
				M3_vals = [ M3_vals opt_out.M3 ];
			else
				M3_vals = [ M3_vals nan];
            end
            fprintf('.')
        end
        disp('.')
        
		results.o1 = opt_out;
		results.e1 = estima;

		figure;
		plot(M1_vals,M3_vals)
		axis([min(M1_vals) max(M1_vals) 0 max(M3_vals)+0.3 ])
		xlabel('M1'); ylabel('M3');


		results.t1.sys = acc_ad;
		results.t1.M1_arr = M1_vals;
		results.t1.M3_arr = M3_vals;
		results.t1.L = L1;
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Iterate through Similar Values when eta_v > eta_w %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if any(attempt_these_tests == 2)

		acc_ad.eta_w = 0.05; acc_ad.eta_v = 0.1;

		M3_vals = [];

		for M1 = M1_vals

			[opt_out,estima] = acc_ad.free_rec_design_pb('Min_M3',M1,M1+2,L1,'verbosity',0);

			if opt_out.problem == 0
				M3_vals = [ M3_vals opt_out.M3 ];
			else
				M3_vals = [ M3_vals nan];
            end
            fprintf('.')
        end
        disp('.')
        
		figure;
		plot(M1_vals,M3_vals)
		axis([min(M1_vals) max(M1_vals) 0 max(M3_vals)+0.3 ])
		xlabel('M1'); ylabel('M3');

		results.t2.sys = acc_ad;
		results.t2.M1_arr = M1_vals;
		results.t2.M3_arr = M3_vals;
		results.t2.L = L1;

	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Iterate through Similar Values when L is "nicer" %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if any(attempt_these_tests == 3)

		L2 = L1;
		L2(end) = 1;

		M3_vals = [];

		for M1 = M1_vals

			[opt_out,estima] = acc_ad.free_rec_design_pb('Min_M3',M1,M1+2,L2,'verbosity',0);

			if opt_out.problem == 0
				M3_vals = [ M3_vals opt_out.M3 ];
			else
				M3_vals = [ M3_vals nan];
            end
            fprintf('.')
        end
        disp('.')

		figure;
		plot(M1_vals,M3_vals)
		axis([min(M1_vals) max(M1_vals) 0 max(M3_vals)+0.3 ])
		xlabel('M1'); ylabel('M3');

		results.t3.sys = acc_ad;
		results.t3.M1_arr = M1_vals;
		results.t3.M3_arr = M3_vals;
		results.t3.L = L2;

	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Iterate through Similar Values when L is ideal %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if any(attempt_these_tests == 4)

		L3 = ones(size(L1));

		M3_vals = [];

		for M1 = M1_vals

			[opt_out,estima] = acc_ad.free_rec_design_pb('Min_M3',M1,M1+2,L3,'verbosity',0);

			if opt_out.problem == 0
				M3_vals = [ M3_vals opt_out.M3 ];
			else
				M3_vals = [ M3_vals nan];
            end
            fprintf('.')
        end
        disp('.')

		figure;
		plot(M1_vals,M3_vals)
		axis([min(M1_vals) max(M1_vals) 0 max(M3_vals)+0.3 ])
		xlabel('M1'); ylabel('M3');

		results.t4.sys = acc_ad;
		results.t4.M1_arr = M1_vals;
		results.t4.M3_arr = M3_vals;
		results.t4.L = L3;

	end

end