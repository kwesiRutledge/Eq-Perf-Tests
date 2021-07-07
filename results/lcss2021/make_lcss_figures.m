% make_lcss_figures.m
%Description:
%	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choosing Figures to Create %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figures_to_create = {'A'}

%%%%%%%%%%%%%%
%% Figure A %%
%%%%%%%%%%%%%%

if any(strcmp(figures_to_create,'A'))
	draw_consistency_set_example1
end