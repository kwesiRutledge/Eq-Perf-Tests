%create_drone_figures.m
%Description:
%   Plots some of the consistency sets for the drone example.

clear all;
close all;
clc;

%% Add Libraries %%

addpath(genpath('../../functions/'))

%% Constants %%

TimeHorizon = 5;
[ lcsas0 , x0 , TimeHorizon , P_target ] = get_differently_loaded_drone_lcsas('TimeHorizon',TimeHorizon);

%% Consistency Sets Plotting %%
L = lcsas0.L;
T = length(L.words{1});
[ n_x , n_u , n_y , n_w , n_v ] = lcsas0.Dimensions();

settings.subset_search_strategy = 'DescendingCardinality';

[ unchecked_knowledge_sequences , sequence_construction_history ] = L.create_belief_sequences_of_length(TimeHorizon);
synthesis_info.UncheckedBeliefSequences = unchecked_knowledge_sequences;

num_knowl_sequences = size(unchecked_knowledge_sequences,2);

[ possible_subsets_of_paths , ~ , choices_as_binary_flags ] = lcsas0.get_feasible_combinations_of_beliefs( unchecked_knowledge_sequences , ...
                                                                                                            'verbosity' , 1 , ... 
                                                                                                            'SkipOptimizationBasedPruning', true );
synthesis_info.PossibleSubsetsOfPaths = possible_subsets_of_paths;

[ possible_subsets_of_paths , choices_as_binary_flags ] = organize_subsets_of_paths( unchecked_knowledge_sequences , choices_as_binary_flags , settings.subset_search_strategy );

%% Get Final Subset

M_Final = possible_subsets_of_paths{1};
[T,num_sequences] = size(M_Final);

M_F_asC = ExternalBehaviorSet.empty;
for M_Final_index = 1:size(M_Final,2)
   %
   M_F_asC(M_Final_index) = ExternalBehaviorSet(lcsas0, M_Final(:,M_Final_index));
   M_F_asPolyhedron(M_Final_index) = M_F_asC(M_Final_index).ToPolyhedron();
end

% Plot at each time
for t = 0:T-1
    figure;
    for word_index = 1:num_sequences
        subplot(1,num_sequences,word_index)
        plot(M_F_asPolyhedron(word_index).projection([1:n_x]+n_x*t))
    end
end