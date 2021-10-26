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

if any(strcmp(figures_to_create,'C'))
	create_cbc_for_toy1
end

if any(strcmp(figures_to_create,'D'))
	create_cbc_for_toy2
end