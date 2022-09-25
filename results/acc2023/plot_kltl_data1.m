%plot_kltl_data1.m
%Description:
%	In this script, we will load some of the data collected from the Crazyflie 2.1 drone and then
%	we will plot it for use in the paper.

%% Constants

mode1_datafilename = 'data/kltl-reach-data1/data/results-mod_pos-measurements-08022022_16_08_18-mode0.txt';

%% Algorithm

% Open the File and Retrieve the data
mode1_file = fopen(mode1_datafilename);
mode1_data = textscan(mode1_file,'%f %f %f %f', 'Delimiter',' ') % Retrieves the four columns and parses them as floats.


% Plot data
figure;
for data_column = 1:length(mode1_data)
	subplot(length(mode1_data),1,data_column)
	plot(mode1_data{data_column})
end


% close all files
fclose(mode1_file);