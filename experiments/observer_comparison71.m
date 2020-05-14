function [results] = observer_comparison71( varargin )
	%observer_comparison71.m
	%Description:
	%	Testing various systems with the algorithms that have been developed so far.
	%

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	% if nargin >= 1
	% 	load_data_flag = varargin{1};
	% end

	% if nargin >= 3
	% 	c_sq.dim_x = varargin{2};
	% 	c_sq.dim_y = varargin{3};
	% end

	% if nargin >= 4
	% 	verbosity = varargin{4};
	% end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	test_name = 'observer_comparison71';
	save_file_name = [ 'results/' test_name '_results.mat'];

	dim = 2;

	L1 = Language([1,2,1,2],[3,4,3,4],[1,1,1,1]);
	T = 4;

	%Create a simple Language Constrainted Switching System
	A1 = [0,1;0.1,-0.05];
	B1 = [0;1];
	C1 = [1,0];
	
	n_x = dim;
	n_u = size(B1,2);
	n_y = size(C1,1);

	eta_v = 0.1; eta_w = 0.2;
	Pv1 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));
	Pw1 = Polyhedron('lb',-eta_w*ones(1,n_x) ,'ub',eta_w*ones(1,n_x));
	% eta_v = 0.3;
	% Pv2 = Polyhedron('lb',-eta_v*ones(1,dim) ,'ub',eta_v*ones(1,dim));

	eta_u = 0; eta_x0 = 0.3;
	P_u = Polyhedron('lb',-eta_u*ones(1,n_u) ,'ub',eta_u*ones(1,n_u));
	P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));

	f1 = eta_w*[0;1];
	f2 = eta_w*[1;0];
	f3 = -f1;
	f4 = -f2;

	aff_dyn_list = [	Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f3,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f4,C1,Pw1,Pv1) ];

	lcsas0 = LCSAS( aff_dyn_list , L1 );

	results.lcsas = lcsas0;

	%%%%%%%%%%%
	%% Tests %%
	%%%%%%%%%%%

	disp(['Beginning ' test_name '.' ])

	%% Construct the Phi Sets for Each Language in Their Entirety
	tic;
	t = 3;
	t = min(t,length(L1.words{1})-1);
	[~, Phi1] = lcsas0.consistent_set(t,Language(L1.words{1}),P_u,P_x0,'fb_method','state','use_proj',false);
	[~, Phi2] = lcsas0.consistent_set(t,Language(L1.words{2}),P_u,P_x0,'fb_method','state','use_proj',false);
	[~, Phi3] = lcsas0.consistent_set(t,Language(L1.words{3}),P_u,P_x0,'fb_method','state','use_proj',false);

	select_y = [ eye(n_x*(t+1)+n_u*t),zeros( n_x*(t+1)+n_u*t , size(Phi1.A,2) - (n_x*(t+1)+n_u*t) ) ];

	results.timing_info.phi_set_construction = toc;
	results.Phi1 = Phi1;
	results.Phi2 = Phi2;
	results.Phi3 = Phi3;

	Phi_list = [Phi1,Phi2,Phi3];

	L1_power = L1.powerset();
	ops0 = sdpsettings('verbose',0);

	%% Compute Initial Tree Nodes %%
	feas_list = [];
	for lang_idx = 1:length(L1_power)
		temp_lang = L1_power(lang_idx);

		%Use Optimization to check if there exists a point that can be explained
		%by all of the words in temp_lang.

		phi = sdpvar(size(Phi1.A,2),temp_lang.cardinality(),'full');

		%Each phi must be contained within each Phi set.
		feas_constr = [];
		for word_idx = 1:temp_lang.cardinality()
			[~,word2phi_idx] = L1.contains( temp_lang.words{word_idx} );
			feas_constr = feas_constr + [ Phi_list(word2phi_idx).A * phi(:,word_idx) <= Phi_list(word2phi_idx).b ] ...
									  + [ Phi_list(word2phi_idx).Ae * phi(:,word_idx) == Phi_list(word2phi_idx).be ];
		end

		%All Phi must have the same measurement.
		matching_obsv_constr = [];
		for word_idx = 2:temp_lang.cardinality()
			matching_obsv_constr = matching_obsv_constr + [ select_y*phi(:,word_idx-1) == select_y*phi(:,word_idx) ];
		end

		diagnostics = optimize(feas_constr+matching_obsv_constr,[],ops0);
		feas_list = [ feas_list , diagnostics ];

		%value(phi)
	end

	results.feas_list = feas_list;

	%%%%%%%%%%%%%%%%%%
	%% Experiment 2 %%
	%%%%%%%%%%%%%%%%%%

	disp('============')
	disp('Experiment 2: Identifying if the range of a linear operator')
	disp('				(on the space of positive vectors) has a range that is all of Rn.')
	disp(' ')

	tic;

	H_y = [Phi1.A; Phi1.Ae;-Phi1.Ae];
	h_y = [Phi1.b; Phi1.be;-Phi1.be];

	H_x = [Phi3.A; Phi3.Ae;-Phi3.Ae];
	h_x = [Phi3.b; Phi3.be;-Phi3.be];

	results.exp2.H_y = H_y;
	results.exp2.h_y = h_y;

	disp('Check to see if H_y'' is nonsingular')
	disp(['rank(H_y'') = ' num2str(rank(H_y')) ]);
	disp(['n_x*(t+1)+n_u*t = ' num2str(n_x*(t+1)+n_u*t)])
	disp(' ')
	disp(['rank(H_y'') == size(H_y'',1) = ',num2str( rank(H_y') == size(H_y',1) ) ])

	H_yp = H_y';
	H_yp_subset = H_yp([1:n_x*(t+1)+n_u*t],:);
	H_yp_rest = H_yp([n_x*(t+1)+n_u*t+1:end],:);

	disp(rank(H_yp_subset))
	disp(['rank(H_yp_subset) == size(H_yp_subset,1) = ',num2str( rank(H_yp_subset) == size(H_yp_subset,1) ) ])

	disp(' ')
	disp('Testing Nonnegative Least Squares')

	A0 = eye(3);
	y0 = [ 3 ; 0 ; -2 ];
	disp(['A = I, y0 = ' num2str(y0') ])
	[x0,fval0,exitflag0,output0] = nonnegative_ls( A0 , y0 );
	disp(['NNLS answer = ' num2str(x0') ])

	A1 = A0;
	y1 = y0; y1(2) = -1;
	disp(['A = I, y1 = ' num2str(y1') ])
	[x1,fval1,exitflag1,output1] = nonnegative_ls( A1 , y1 );
	disp(['NNLS answer = ' num2str(x1') ])	

	%Verify that the full state space is recoverable with this matrix using nonnegative numbers.
	disp('Verifying that the range of the operator H_y'' is the full state space.')
	col_count = 0;
	for col_idx = 1:size(H_yp_subset,2)
		col_out = H_yp_subset(:,col_idx);
		temp_H = H_yp_subset(:,[1:size(H_yp_subset,2)] ~= col_idx);	

		[xi,fvali,exitflagi,outputi] = nonnegative_ls( temp_H , -col_out );
		fval_list(col_idx) = fvali;
		if fvali < 1e-4
			col_count = col_count + 1;
		end
	end
	disp(['- Are all columns capable of constructing each other'])
	disp(['using POSITIVE combination coefficients? ' num2str(col_count == size(H_yp,2)) ])
	results.exp2.col_count = col_count;
	results.exp2.H_yp_subset = H_yp_subset;
	results.exp2.fval_list = fval_list;

	results.range_checking_time = toc;

	%Verify that the full state space is recoverable with this matrix using nonnegative numbers.
	col_count = 0;
	for col_idx = 1:size(H_yp_subset,2)
		col_out = H_yp_subset(:,col_idx);
		temp_H = H_yp_subset(:,[1:size(H_yp_subset,2)] ~= col_idx);	

		[xi,fvali,exitflagi,outputi] = constrained_nonnegative_ls( temp_H , -col_out , ...
																	[] , [] , ...
																	H_yp_rest , zeros(size(H_yp_rest,1),1) , ...
																	'verbosity', 0 );
		fval_list(col_idx) = fvali;
		if fvali < 1e-4
			col_count = col_count + 1;
		end
	end
	disp(['- Are all columns capable of constructing each other'])
	disp(['using POSITIVE AND CONSTRAINED combination coefficients? ' num2str(col_count == size(H_yp,2)) ])
	
	results.exp2.constrained_fval_list = fval_list;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Save The Larger Variables to a Data File %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%save([ save_file_name '.mat'])

end