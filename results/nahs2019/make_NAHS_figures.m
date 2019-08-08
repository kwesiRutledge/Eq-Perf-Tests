%make_NAHS_figures.m

cd ../..

%Run the Comparision of Worst-CAse Language vs. Prefix-Based Language
ot1 = observer_tests(28);

%Run the Zonotope Example
% zonotope_template = observer_tests(49);

%Run the simple lane keeping controller example
ot2 = observer_tests(37);

%Run the Multiple Agent Control Example
multi_agent = observer_tests([50,50,50],{'2,2','3,3','4,4'});

%