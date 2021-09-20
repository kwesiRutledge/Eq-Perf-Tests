%make_NAHS_like_video.m

clear all;
close all;
clc;

%Get Data
if exist('NAHS_video_data0.mat') ~= 2 
    cd ../..
    %If the NAHS_video trajectory does not exist, then run the simulation
    %again.
    multi_agent_data = observer_tests([50],{'3,3,5,true'});
    cd results/nahs2019
    save('NAHS_video_data0.mat','multi_agent_data')
else
    load('NAHS_video_data0.mat')
end

%Process this data into useable trajectory data for OpenRave
quadrotor_xyzQuat_trajectories = {};
num_quadrotors = multi_agent_data{1}.params.cube_x*multi_agent_data{1}.params.cube_y;
for quadrotor_idx = 1:num_quadrotors
    %Quadrotor Postitions
    quadrotor_i_xpos = multi_agent_data{1}.x_trajectory(quadrotor_idx,:);
    quadrotor_i_ypos = multi_agent_data{1}.x_trajectory(num_quadrotors+quadrotor_idx,:);
    
    switch mod(quadrotor_idx,2)
        case 0
            quadrotor_i_zpos = -1.0*ones(size(quadrotor_i_xpos));
        case 1
            quadrotor_i_zpos = 1.0*ones(size(quadrotor_i_xpos));
        otherwise
            quadrotor_i_zpos = NaN(size(quadrotor_i_xpos));
    end
    
    %Save to new variable
    quadrotor_xyzQuat_trajectories{quadrotor_idx} = [ quadrotor_i_xpos', quadrotor_i_ypos', quadrotor_i_zpos' , ones(length(quadrotor_i_xpos),1), zeros(length(quadrotor_i_xpos),3) ];
    
end

save('NAHS_video_data0_forOpenRave.mat','quadrotor_xyzQuat_trajectories')



