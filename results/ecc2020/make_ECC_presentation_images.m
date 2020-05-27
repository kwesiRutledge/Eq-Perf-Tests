%make_ECC_presentation_images.m
%Description:
%	This function is meant to create the images that are used in the ECC presentation.

cd ../..

clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Bare Graphs Using the Saved Gains %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% results = observer_tests(60,{[	' ''c_sq'', 2,1 ,' ...
% 			 					' ''L'' ,Language([1,2,2,1,1,1],[1,2,1,1,1,1],[1,1,2,1,1,1]),' ...
% 			 					' ''disturb_bounds'',0.35,0.2,' ...
% 			 					' ''eta_u'', 50 ,' ...
% 			 					' ''load_data_flag'', true, ' ...
% 			 					' ''save_file_name'', ''data/ecc/oc60_interm_results_2x1drones_23_57_18'' ']})

% results = observer_tests(74,{[	' ''c_sq'', 2,1 ,' ...
% 			 					' ''L'' ,Language([1,2,2,1,1,1],[1,2,1,1,1,1],[1,1,2,1,1,1]),' ...
% 			 					' ''disturb_bounds'',0.35,0.2,' ...
% 			 					' ''eta_u'', 50 ,' ...
% 			 					' ''load_data_flag'', true, ' ...
% 			 					' ''save_file_name'', ''data/ecc/oc60_interm_results_2x1drones_23_57_18'' ']})

% results = observer_tests(74,{[	' ''c_sq'', 2,1 ,' ...
% 								' ''L'' ,Language([1,2,2,1,1],[1,2,1,1,1],[1,1,2,1,1]),' ...
% 								' ''disturb_bounds'',0.35,0.2,' ...
% 								' ''eta_u'', 50 ,' ...
% 								' ''load_data_flag'', true, ' ...
% 			 					' ''save_file_name'', ''data/ecc/oc60_interm_results_2x1drones_23_09_28'' ']})

results = observer_tests(74,{[	' ''c_sq'', 2,1 ,' ...
								' ''L'' ,Language([1,2,2,1],[1,2,1,1],[1,1,2,1]),' ...
								' ''disturb_bounds'',0.35,0.2,' ...
								' ''eta_u'', 50 ,' ...
								' ''load_data_flag'', true, ' ...
 			 					' ''save_file_name'', ''data/ecc/oc74_interm_results_2x1drones_19_41_24'' ']})