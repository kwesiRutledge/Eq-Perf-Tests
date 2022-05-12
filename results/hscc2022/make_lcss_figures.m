% make_lcss_figures.m
%Description:
%	Creates the figures which were used in the submission to HSCC 2022.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choosing Figures to Create %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figures_to_create = {'I'}

%Figure Options:
%   A = Draw Consistency Sets (overlapping) for the opposing rotation system
%   B = Create summary figure.
%   C = Creates Controller for the Opposing Rotation System ?
%   D = Creates Controller and Plots Closed Loop Response for the Similar Rotation System
%   E = Creates figure showing the performance of the Differently Loaded Drone Performance
%   F = Creates alternative plot for similar rotation system.
%   G = Draws consistency sets for a specific instant in time for opposing rotation.
%   H = Draws a single sequence of states and labels the hypothesis on each one.


%% Iterate through elements of figures_to_create

for figure_index = 1:length(figures_to_create)
    
    switch figures_to_create{figure_index}
        %Figure A
        case 'A'
            draw_consistency_set_example1
        case 'B'
            draw_summary_figure
        case 'C'
            create_cbc_for_toy1
        case 'D'
            create_cbc_for_toy2
        case 'D-Simple'

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Create Figure D: Plot of Similar Rotation Performance WITHOUT finding Controller %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
            % Load good data for this plot. If this data does not
            % exist, then you need to run create_cbc_for_drone first
            % and then copy the name of the data file it creates below.
            load('data/toy2_data_05Oct2021-1833')

            % Constants
            colorCode = {'cyan','magenta'};
            lineWidth = 2;

            target_color = 'yellow';
            target_alpha = 0.5;

            fs = 20

            figure('DefaultAxesFontSize',fs);

            hold on;
            plot(lcsas0.X0)
            plot(P_target,'color',target_color,'alpha',target_alpha)

            for simulation_index = 1:15
                [ x_0_t, u_0_tm1 , y_0_t , sig ] = toy2_controller.simulate_1run();

                % for t = 0:TimeHorizon
                %   x_t = x_0_t(:,t+1);
                %   scatter(x_t(1),x_t(2))
                % end
                plot(x_0_t(1,:),x_0_t(2,:),'Color',colorCode{sig(1)},'LineWidth',lineWidth )
                

                %Save data
                results.SimulationData = [ results.SimulationData ; struct('x_0_t',x_0_t,'u_0_tm1',u_0_tm1) ];

            end
            axis([-0.5,12.5,-6,6.5])

            xlabel('$x_1$','Interpreter','latex')
            ylabel('$x_2$','Interpreter','latex')

            saveas(gcf,'images/toy2_runs','epsc')
            saveas(gcf,'images/toy2_runs','png')

        case 'E'
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Create Figure E: Plot of Drone Performance %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
            % Load good data for this plot. If this data does not
            % exist, then you need to run create_cbc_for_drone first
            % and then copy the name of the data file it creates below.
            load('data/drone_data_30Oct2021-0454')

            % Constants
            colorCode = {'cyan','magenta'};
            lineWidth = 2;

            target_color = 'yellow';
            target_alpha = 0.5;

            fs = 20;

            figure('DefaultAxesFontSize',fs);

            hold on;
            plot(lcsas0.X0)
            plot(P_target,'color',target_color,'alpha',target_alpha)

            for simulation_index = 1:15
                [ x_0_t, u_0_tm1 , y_0_t , sig ] = drone_controller.simulate_1run();


                plot(x_0_t(1,:),x_0_t(2,:),'Color',colorCode{sig(1)},'LineWidth',lineWidth)

                %Save data
                results.SimulationData = [ results.SimulationData ; struct('x_0_t',x_0_t,'u_0_tm1',u_0_tm1) ];

            end
            axis([0.5,5,-4,10])

            xlabel('$x_1$','Interpreter','latex')
            ylabel('$x_2$','Interpreter','latex')

            saveas(gcf,'images/drone_runs2','epsc')
            saveas(gcf,'images/drone_runs2','png')

        case 'F'
            % Create alternative plof for toy2
            
            % Load good data for this plot. If this data does not
            % exist, then you need to run create_cbc_for_drone first
            % and then copy the name of the data file it creates below.
            load('data/toy2_data_05Oct2021-1833')
            
            Dim = lcsas0.X0.Dim;
            T = TimeHorizon;
            
            n_simulations = 15;

            colorCode = {'cyan','magenta'};
            lineWidth = 3;

            % Create data
            for simulation_index = 1:n_simulations
                [ x_0_t, u_0_tm1 , y_0_t , sig ] = toy2_controller.simulate_1run();
                x_0_t_arr{simulation_index} = x_0_t;
                sig_arr{simulation_index} = sig;
            end

            axesOptions = {[-0.5,TimeHorizon+0.5,-0.5,12],[-0.5,TimeHorizon+0.5,-5,8]}
            % Create figure
            figure;

            for dim_index = 1:Dim
                subplot(Dim,1,dim_index)
                
                hold on;
                plot(0,lcsas0.X0.V(dim_index))
                plot(Polyhedron('lb',T,'ub',T)*P_target.projection(dim_index))

                % for t = 0:TimeHorizon
                % 	x_t = x_0_t(:,t+1);
                % 	scatter(x_t(1),x_t(2))
                % end
                for simulation_index = 1:n_simulations
                    sig_for_si = sig_arr{simulation_index}(1);
                    plot([0:T],x_0_t_arr{simulation_index}(dim_index,:),'Color',colorCode{sig_for_si},'LineWidth',lineWidth)
                end

                axis(axesOptions{dim_index})
            end
            %axis([-0.5,12.5,-6,6.5])

            %Save data
            results.SimulationData = [ results.SimulationData ; struct('x_0_t',x_0_t,'u_0_tm1',u_0_tm1) ];


            saveas(gcf,'images/toy2_runs_v2','epsc')
            saveas(gcf,'images/toy2_runs_v2','png')

        case 'G'
            % Create a plot showing the difference between the multiple consistency sets for the similar
            % rotation system.
            draw_consistency_set_example3
        case 'H'
            %Plot a sequence of states
            load('data/toy2_data_21Dec2021-1607')

            % Simulate the system ONCE
            [ x_0_t, u_0_tm1 , y_0_t , sig ] = toy2_controller.simulate_1run();
            
            extreme_vector = [ min([x_0_t(1,:),P_target.V(:,1)']) , max([x_0_t(1,:),P_target.V(:,1)']) , min([x_0_t(2,:),P_target.V(:,2)']) , max([x_0_t(2,:),P_target.V(:,2)']) ];
            marker_size = 72; %Default is 36

            figure;
            hold on;
            plot(P_target)
            plot(x_0_t(1,:),x_0_t(2,:), 'LineWidth', 2.0)
            % for t = 0:TimeHorizon
            %     % Label each point with their time
            %     text(x_0_t(1,t+1),x_0_t(2,t+1) - 0.2,['$$\uparrow t = ' num2str(t) '$$'] , 'Interpreter' , 'latex' )
            % end

            for t = 0:TimeHorizon-1
                % Label each point with their time
                %text(x_0_t(1,t+1),x_0_t(2,t+1) - 0.2,['$$\uparrow t = ' num2str(t) '$$'] , 'Interpreter' , 'latex' )


                %Label each point with the hypothesis
                mu_t = toy2_controller.b_hist(t+1);
                mu_string = ['$$\mu_{' num2str(t) '} = \{' ];
                for mu_index = 1:mu_t.cardinality()
                    word_i = mu_t.words{mu_index};

                    %mu_string = [ mu_string '\mathbf{m}^{(' num2str(word_i(1)) ')}' ];
                    mu_string = [ mu_string num2str(word_i(1)) ];

                    % Add comma if this is not the last word_i in mu_t
                    if mu_index ~= mu_t.cardinality()
                        mu_string = [ mu_string ',' ];
                    end
                end
                %Finish latex string
                mu_string = [ mu_string ' \}$$' ];

                text(x_0_t(1,t+1),x_0_t(2,t+1) - 0.2 , mu_string , 'Interpreter' , 'latex' )
            end
            
            %plot(x_0_t(1,:),x_0_t(2,:))

            axis(extreme_vector + [ -0.5 , 0.5 , -0.5 , 0.5 ])

            saveas(gcf,'images/toy2_detection_diagram','epsc')
            saveas(gcf,'images/toy2_detection_diagram','png')


        case 'I'
            draw_consistency_set_drone_example

        end
            
            
end

