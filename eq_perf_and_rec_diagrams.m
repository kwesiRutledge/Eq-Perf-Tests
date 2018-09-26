% eq_perf_and_rec_diagrams.m

clear all close all;

%% Create Equalized Performance Diagram

xh_t = [1;2];
x_t = xh_t + unifrnd(-1,1,[2,1]);

M = 1;

fs = 20 ; %Font Size for Plots

figure;
subplot(1,2,1)
hold on;
scatter(x_t(1),x_t(2),'ro')
scatter(xh_t(1),xh_t(2),'bx')
rectangle('Position',[ xh_t' - ones(size(xh_t))'*M 2*M 2*M ], ...
            'EdgeColor','b',...
            'LineWidth',3)

axis([-1 4 -1 4])
grid on
xlabel('$x_1(t)$','Interpreter','latex','FontSize',20)
ylabel('$x_2(t)$','Interpreter','latex','FontSize',20)

% Create time t+1
xh_tp1 = [2;1];
x_tp1 = xh_tp1 + unifrnd(-1,1,[2,1]);

subplot(1,2,2)
hold on;
scatter(x_tp1(1),x_tp1(2),'ro')
scatter(xh_tp1(1),xh_tp1(2),'bx')
rectangle('Position',[ xh_tp1' - ones(size(xh_tp1))'*M 2*M 2*M ], ...
            'EdgeColor','b',...
            'LineWidth',3)

axis([-1 4 -1 4])
grid on
xlabel('$x_1(t+1)$','Interpreter','latex','FontSize',fs)
ylabel('$x_2(t+1)$','Interpreter','latex','FontSize',fs)


%% Equalized Recovery Diagram

%Final Step Enforces guarantee.
xh_tpT = xh_tp1;
x_tpT = x_tp1;

%Intermediate Steps have a different guarantee.
M2 = 1.5;

xh_tpk = xh_t*0.5 + xh_tpT*0.5;
x_tpk = xh_tpk + unifrnd(-M2,M2,2,1);

figure;

subplot(1,3,1)
hold on;
scatter(x_t(1),x_t(2),'ro')
scatter(xh_t(1),xh_t(2),'bx')
rectangle('Position',[ xh_t' - ones(size(xh_t))'*M 2*M 2*M ], ...
            'EdgeColor','b',...
            'LineWidth',3)

axis([-1 4 -1 4])
grid on
xlabel('$x_1(t)$','Interpreter','latex','FontSize',fs)
ylabel('$x_2(t)$','Interpreter','latex','FontSize',fs)

subplot(1,3,2)
hold on;
scatter(x_tpk(1),x_tpk(2),'ro')
scatter(xh_tpk(1),xh_tpk(2),'bx')
rectangle('Position',[ xh_tpk' - ones(size(xh_tp1))'*M2 2*M2 2*M2 ], ...
            'EdgeColor','b',...
            'LineWidth',3)

axis([-1 4 -1 4])
grid on
xlabel('$x_1(t+k)$','Interpreter','latex','FontSize',fs)
ylabel('$x_2(t+k)$','Interpreter','latex','FontSize',fs)

subplot(1,3,3)
hold on;
scatter(x_tpT(1),x_tpT(2),'ro')
scatter(xh_tpT(1),xh_tpT(2),'bx')
rectangle('Position',[ xh_tpT' - ones(size(xh_tp1))'*M 2*M 2*M ], ...
            'EdgeColor','b',...
            'LineWidth',3)

axis([-1 4 -1 4])
grid on
xlabel('$x_1(t+T)$','Interpreter','latex','FontSize',fs)
ylabel('$x_2(t+T)$','Interpreter','latex','FontSize',fs)

