% function tests = testBeliefGraph
% 	%disp(localfunctions)
% 	tests = functiontests(localfunctions);

% function test1_BeliefGraph(testCase)
% 	%Description:
% 	%	Testing that the two BeliefGraph objects from buildTwoGraphs have the proper relationship.
% 	%	(1 is covered by 2, and 2 is NOT covered by 1)

% 	[~,~,~,bg1,bg2] = buildTwoGraphs();

% 	assert( bg1.is_covered_by(bg2) & (~bg2.is_covered_by(bg1)) )

% function test2_BeliefGraph(testCase)
% 	%Description:
% 	%	Testing that the two BeliefGraph objects from buildTwoGraphs contain themselves.

% 	[~,~,~,bg1,bg2] = buildTwoGraphs();

% 	assert( bg1.is_covered_by(bg1) & bg1.is_covered_by(bg1) )

% function test3_BeliefGraph(testCase)
% 	%Description:
% 	%	Testing the is_abstracted_by() function.
% 	%	For first pair of BeliefGraphs from buildTwoGraphs(), it should be the case that they both
% 	%	abstract one another.

% 	[~,~,~,bg1,bg2] = buildTwoGraphs(1);	

% 	assert( bg2.is_abstracted_by(bg1) & bg1.is_abstracted_by(bg2) )

% function test4_BeliefGraph(testCase)
% 	%Description:
% 	%	Testing the _is_abstracted_by() function.
% 	%	For the SECOND pair of BeliefGraphs from buildTwoGraphs(), it should be the case
% 	%	that bg2 is abstracted by bg1 but NOT the other way around.

% 	[~,~,~,bg1,bg2] = buildTwoGraphs(2);	

% 	assert( bg2.is_abstracted_by(bg1) & (~bg1.is_abstracted_by(bg2)) )	

% function [L1,lcsas0,empty_bg,bg1,bg2] = buildTwoGraphs(varargin)

% 	%% Input Processing %%

% 	if nargin >= 1
% 		version_num = varargin{1};
% 	end

% 	if ~exist('version_num')
% 		version_num = 1;
% 	end



% 	%%%%%%%%%%%%%%%
% 	%% Algorithm %%
% 	%%%%%%%%%%%%%%%

% 	switch version_num
% 		case 1
% 			L1 = Language([1,2,1,2],[3,4,3,4],[5,1,1,1]);
% 			T = 4;

% 			%Create a simple Language Constrainted Switching System
% 			A1 = [0,1;0.1,-0.05];
% 			B1 = [0;1];
% 			C1 = [1,0];

% 			n_x = size(A1,1);
% 			n_u = size(B1,2);
% 			n_y = size(C1,1);

% 			eta_v = 0.1; eta_w = 0.2;
% 			Pv1 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));
% 			Pw1 = Polyhedron('lb',-eta_w*ones(1,n_x) ,'ub',eta_w*ones(1,n_x));
% 			eta_v = 0.3;
% 			Pv2 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));

% 			eta_u = 0; eta_x0 = 0.3;
% 			P_u = Polyhedron('lb',-eta_u*ones(1,n_u) ,'ub',eta_u*ones(1,n_u));
% 			P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));

% 			f1 = 0*[0;1];
% 			f2 = 0*[1;0];
% 			f3 = -f1;
% 			f4 = -f2;

% 			aff_dyn_list = [	Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1), ...
% 								Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1), ...
% 								Aff_Dyn(A1,B1,f3,C1,Pw1,Pv1), ...
% 								Aff_Dyn(A1,B1,f4,C1,Pw1,Pv1), ...
% 								Aff_Dyn(A1,B1,f1,C1,Pw1,Pv2) ];

% 			lcsas0 = LCSAS( aff_dyn_list , L1 , 'X0' , P_x0 );
% 			lcsas0.X0 = P_x0;

% 			empty_bg = BeliefGraph( lcsas0 , P_u , P_x0 , 'return_empty' , true );

% 			%Hand construct BG1
% 			bg1 = empty_bg;
% 			for tau = 0:T-1
% 				temp_bn = BeliefNode( L1 , tau );
% 				bg1.N = [ bg1.N , temp_bn ];
% 			end

% 			for tau = 1:T-1
% 				temp_edge = [tau,tau+1];
% 				bg1.E = [bg1.E; temp_edge];
% 			end

% 			bg2 = BeliefGraph( lcsas0 , P_u , P_x0 , ...
% 								'verbosity', 0 , 'accel_flag', true , ...
% 								'use_proj_flag', true , 'use_unobservability_checks', true );
% 		case 2

% 			L1 = Language([1,2,1,2],[3,4,3,4],[5,1,1,1]);
% 			T = 4;

% 			%Create a simple Language Constrainted Switching System
% 			A1 = [0,1;0.1,-0.05];
% 			B1 = [0;1];
% 			C1 = [1,0];

% 			n_x = size(A1,1);
% 			n_u = size(B1,2);
% 			n_y = size(C1,1);

% 			eta_v = 0.2; eta_w = 0.2;
% 			Pv1 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));
% 			Pw1 = Polyhedron('lb',-eta_w*ones(1,n_x) ,'ub',eta_w*ones(1,n_x));
% 			eta_v = 0.3;
% 			Pv2 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));

% 			eta_u = 0; eta_x0 = 0.3;
% 			P_u = Polyhedron('lb',-eta_u*ones(1,n_u) ,'ub',eta_u*ones(1,n_u));
% 			P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));

% 			f1 = eta_w*[0;1];
% 			f2 = eta_w*[1;0];
% 			f3 = -f1;
% 			f4 = -f2;

% 			aff_dyn_list = [	Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1), ...
% 								Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1), ...
% 								Aff_Dyn(A1,B1,f3,C1,Pw1,Pv1), ...
% 								Aff_Dyn(A1,B1,f4,C1,Pw1,Pv1), ...
% 								Aff_Dyn(A1,B1,f1,C1,Pw1,Pv2) ];

% 			lcsas0 = LCSAS( aff_dyn_list , L1 , 'X0' , P_x0 );
% 			lcsas0.X0 = P_x0;

% 			empty_bg = BeliefGraph( lcsas0 , P_u , P_x0 , 'return_empty' , true );

% 			%Hand construct BG1
% 			bg1 = empty_bg;
% 			for tau = 0:T-1
% 				temp_bn = BeliefNode( L1 , tau );
% 				bg1.N = [ bg1.N , temp_bn ];
% 			end

% 			for tau = 1:T-1
% 				temp_edge = [tau,tau+1];
% 				bg1.E = [bg1.E; temp_edge];
% 			end

% 			bg2 = BeliefGraph( lcsas0 , P_u , P_x0 , ...
% 								'verbosity', 0 , 'accel_flag', true , ...
% 								'use_proj_flag', true , 'use_unobservability_checks', true );
% 		otherwise
% 			error('Unrecognized version_num!')
% 	end

classdef testBeliefGraph < matlab.unittest.TestCase
 
    properties
        BeliefGraph1;
        BeliefGraph2;
        BeliefGraph3;
        Lang1;
        System1;
        System3;
    end
 
    methods(TestMethodSetup)
        function createBeliefGraphs(testCase)
        	%FigurePropertiesTest.createBeliefGraphs
        	%Description:
        	%	Setup three Simple BeliefGraphs
        	%	Automatically called when this test is run.

        	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        	%% Define Switching Language %%
        	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        	testCase.Lang1 = Language([1,1,1,1],[2,2,2,2],[3,1,1,1]);
        	T = length(testCase.Lang1.words{1});

        	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        	%% Create Simple 1-Dimensional System (v1) %%
        	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        	%Create a simple Language Constrainted Switching System
			A1 = 0.9;
			B1 = 0;
			C1 = 1;

			B3 = 1;

			n_x = size(A1,1);
			n_u = size(B1,2);
			n_y = size(C1,1);

			eta_v = 0.2; eta_w = 0.2;
			Pv1 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));
			Pw1 = Polyhedron('lb',-eta_w*ones(1,n_x) ,'ub',eta_w*ones(1,n_x));
			eta_v = 0.3;
			Pv2 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));

			eta_u = 0; eta_x0 = 0.3;
			P_u = Polyhedron('lb',-eta_u*ones(1,n_u) ,'ub',eta_u*ones(1,n_u));
			P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));

			f1 = eta_w*0.0;
			f2 = -eta_w*0.0;

			aff_dyn_list1 = [	Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1), ...
								Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1), ...
								Aff_Dyn(A1,B1,f1,C1,Pw1,Pv2) ];

			lcsas_v1 = LCSAS( aff_dyn_list1 , testCase.Lang1 , 'X0' , P_x0 );
			lcsas_v1.X0 = P_x0;

			testCase.System1 = lcsas_v1;

			%Find BeliefGraph
			testCase.BeliefGraph1 = BeliefGraph( testCase.System1 , P_u , P_x0 , ...
													'verbosity', 0 , 'accel_flag', true , ...
													'use_proj_flag', true , 'use_unobservability_checks', true );

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        	%% Create Simple 1-Dimensional System (v2) %%
        	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        	empty_bg = BeliefGraph( testCase.System1 , P_u , P_x0 , 'return_empty' , true );

			%Hand construct BG1
			bg2 = empty_bg;
			for tau = 0:T-1
				temp_bn = BeliefNode( testCase.Lang1 , tau );
				bg2.N = [ bg2.N , temp_bn ];
			end

			for tau = 1:T-1
				temp_edge = [tau,tau+1];
				bg2.E = [bg2.E; temp_edge];
			end

			testCase.BeliefGraph2 = bg2;

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        	%% Create Simple 1-Dimensional System (v3) %%
        	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        	%Modify Switched System
        	eta_w = 0.3;
			Pw1 = Polyhedron('lb',-eta_w*ones(1,n_x) ,'ub',eta_w*ones(1,n_x));

        	f1 = eta_w*1.5;
			f2 = -eta_w*1.5;

			aff_dyn_list3 = [	Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1), ...
								Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1), ...
								Aff_Dyn(A1,B1,f1,C1,Pw1,Pv2) ];

			lcsas_v3 = LCSAS( aff_dyn_list3 , testCase.Lang1 , 'X0' , P_x0 );
			lcsas_v3.X0 = P_x0;
			testCase.System3 = lcsas_v3;

			testCase.BeliefGraph3 = BeliefGraph( testCase.System3 , P_u , P_x0 , ...
													'verbosity', 0 , 'accel_flag', true , ...
													'use_proj_flag', true , 'use_unobservability_checks', true );
        end
    end
 
    methods(TestMethodTeardown)
        % function closeFigure(testCase)
        %     disp('Fin.')
        % end
    end
 
    methods(Test)
 
        function trivial_covering1(testCase)
			%Description:
			%	Testing that the two BeliefGraph objects from buildTwoGraphs have the proper relationship.
			%	(1 is covered by 2, and 2 is NOT covered by 1)

			bg1 = testCase.BeliefGraph1;
			bg2 = testCase.BeliefGraph2;

			assert( bg2.is_covered_by(bg1) & (~bg1.is_covered_by(bg2)) )
		end

		function self_covering1(testCase)
			%Description:
			%	Testing that the first two BeliefGraph objects from buildTwoGraphs contain themselves.

			bg1 = testCase.BeliefGraph1;
			bg2 = testCase.BeliefGraph2;

			assert( bg1.is_covered_by(bg1) & bg1.is_covered_by(bg1) )
		end

		function belief_graph_abstraction1(testCase)
			%Description:
			%	Testing the is_abstracted_by() function.
			%	For first pair of BeliefGraphs from buildTwoGraphs(), it should be the case that they both
			%	abstract one another.

			bg1 = testCase.BeliefGraph1;
			bg2 = testCase.BeliefGraph2;

			assert( bg2.is_abstracted_by(bg1) & bg1.is_abstracted_by(bg2) )
 		end

 		function belief_graph_abstraction2(testCase)
			%Description:
			%	Testing the _is_abstracted_by() function.
			%	For the SECOND pair of BeliefGraphs from buildTwoGraphs(), it should be the case
			%	that bg3 is abstracted by bg2 but NOT the other way around.

			bg2 = testCase.BeliefGraph2;
			bg3 = testCase.BeliefGraph3;

			% figure;
			% subplot(1,2,1)
			% bg2.plot()

			% subplot(1,2,2)
			% bg3.plot()

			% disp( bg3.N(end).subL )
			% bg3.N(end).c_set.V

			assert( bg3.is_abstracted_by(bg2) & (~bg2.is_abstracted_by(bg3)) )	
		end
 
    end
 
end