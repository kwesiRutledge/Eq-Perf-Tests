%make_NAHS_figures.m

cd ../..

%Run the Comparision of Worst-CAse Language vs. Prefix-Based Language
ot1 = observer_tests(28);

%Run the Zonotope Example
% zonotope_template = observer_tests(49);

%Run the simple lane keeping controller example
ot2 = observer_tests(37);

%Run the Multiple Agent Control Example
multi_agent = observer_tests([50,50,50],{'2,2,5,true','3,3,5,true','4,4,5,true'});

for experim_idx = 1:3
	SolverTime(experim_idx,1) = multi_agent{experim_idx}.history.oo3.solvertime;
	YALMIPTime(experim_idx,1) = multi_agent{experim_idx}.history.oo3.yalmiptime;
	CubeX(experim_idx,1) = multi_agent{experim_idx}.params.cube_x;
	CubeY(experim_idx,1) = multi_agent{experim_idx}.params.cube_y;
end

ExperimentNum = [1:3]';
T = table(ExperimentNum,CubeX,CubeY,YALMIPTime,SolverTime)

cd results/nahs2019
save nahs_results4_22aug2019.mat