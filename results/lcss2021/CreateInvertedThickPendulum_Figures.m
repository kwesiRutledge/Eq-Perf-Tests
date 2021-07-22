%CreateInvertedThickPendulum_Figures.m
%Description:
%	Plotting stuff for research update.

clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Include Relevant Libraries %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('PusherSlider') == 0
    %If the class does not exist on the path,
    %then add the systems directory to the path.
    addpath(genpath('~/Documents/Michigan/Research/systemRepository/systems'));
end

%%%%%%%%%%%%%%%
%% Constants %%
%%%%%%%%%%%%%%%

itp1 = InvertedThickPendulum();
itp1.CoMx_rel = -0.75;

itp2 = InvertedThickPendulum();
itp2.CoMx_rel = 0;


FileName1 = 'itp3_show';
% vidObj = VideoWriter(FileName1,'MPEG-4');

%%%%%%%%%%%%%%%%%%%%
%% Create Image 1 %%
%%%%%%%%%%%%%%%%%%%%

ax1 = [-2 2 -2 2];

figure;
subplot(1,2,1)
[all_elts,pend_elts,com_elt,state_elt] = itp1.Show('CoMMarkerSize',5*72,'PivotMarkerSize',2.5*72,'StateMarkerSize',2.5*72)
axis(ax1)
set(gca, 'XAxisLocation', 'origin')
set(gca, 'YAxisLocation', 'origin')
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
box on


% delete(com_elt)

subplot(1,2,2)
[all_elts,pend_elts,com_elt,state_elt] = itp2.Show('CoMMarkerSize',5*72,'PivotMarkerSize',2.5*72,'StateMarkerSize',2.5*72)
axis(ax1)
set(gca, 'XAxisLocation', 'origin')
set(gca, 'YAxisLocation', 'origin')
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
box on

set(gcf,'Units','Normalized','Position',[0,0,1,1])

saveas(gcf,'images/two_options_figure','epsc')
saveas(gcf,'images/two_options_figure','png')
