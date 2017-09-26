% comp_sys.m
% Comparing System Matrices

clear all;
close all;
clc;

%% Tiancheng's System Matrices

% Constants

ki = 0.2890; % already divided by m
f_0 = 0.6558; % already divided by m
f_1 = 0.0058; % already divided by m
T = 1/25; % controller period

% Continuous matrices
A = [-f_1, 0, 0;
     -1, 0, 1;
      0, 0, 0;];
B = [ki;
     0;
     0;];
F = [-f_0;
      0;
      0;];
E = [0;
     0;
     1;];
C = [1, 0, 0;
     0, 1, 0];

% Convert to discret time
A_s = @(s) expm(s*A);
Ad = A_s(T);
Bd = integral(A_s, 0, T, 'ArrayValued', true) * B;
Fd = integral(A_s, 0, T, 'ArrayValued', true) * F;
Ed = integral(A_s, 0, T, 'ArrayValued', true) * E;

% Ad = [Aaa, Aau;
%       Aua, Auu];
Aaa = Ad(1:2,1:2);
Aau = Ad(1:2,3);
Aua = Ad(3,1:2);
Auu = Ad(3,3);

Ba = Bd(1:2,1);
Bu = Bd(3,1);
Fa = Fd(1:2,1);
Fu = Fd(3,1);

Ea = Ed(1:2,1);
Eu = Ed(3,1);

%% My System Matrices

load('data/system_examples/acc_m1_systems.mat')

%% Use Constants from https://github.com/pettni/cps-inv/tree/master/acc_3d_nonconvex

% car dynamics
con.mass = 1370;
con.f0 = 51;
con.f1 = 1.2567;
con.f2 = 0.4342;
con.g = 9.82;

% linearization
con.lin_speed = 10;
con.f0_bar = con.f0 - con.f2*con.lin_speed^2;
con.f1_bar = con.f1 + 2*con.f2*con.lin_speed;
con.dT = 0.5;

% Matrices

s3.A = [ (1/con.mass)*[ -con.f1_bar 0 0 ] ;
		 -1 0 1 ;
		 0 0 0 ];

s3.B = [ (1/con.mass);
		 0 ;
		 0 ];

s3.C = C;

s3.E = E;

s3.K = F;

clear A_s

A_s = @(s) expm(s*s3.A);
s3.Ad = A_s(T);
s3.Bd = integral(A_s, 0, T, 'ArrayValued', true) * s3.B;
s3.Kd = integral(A_s, 0, T, 'ArrayValued', true) * s3.K;
s3.Ed = integral(A_s, 0, T, 'ArrayValued', true) * s3.E;

%% Compare A,B

disp('==========')
disp('A Matrices')
% disp( ...
% [ Ad Inf*ones(3,1) acc_dsys.A Inf*ones(3,1) s3.Ad ] )

method_names = { 'Tiancheng''s' , 'Mine' , 'Petter''s Constants' };
A_cell = { Ad , acc_dsys.A , s3.Ad };
B_cell = { Bd , acc_dsys.B , s3.Bd };
E_cell = { Ed , acc_dsys.E , s3.Ed };
F_cell = { Fd , acc_dsys.F , s3.Kd };

comparing = {'A','B','E','F'};

for mat_num = 1:length(comparing)
	disp('======================')
	disp([comparing{mat_num} ' Matrices'])

	for method_num = 1 : length(method_names)
		disp(method_names{method_num})
		if strcmp(comparing{mat_num},'A') 
			disp(A_cell{method_num})
		end

		if strcmp(comparing{mat_num},'B')
			disp(B_cell{method_num})
		end

		if strcmp(comparing{mat_num},'E')
			disp(E_cell{method_num})
		end

		if strcmp(comparing{mat_num},'F')
			disp(F_cell{method_num})
		end

	end

end

% Save Some of these matrices
acc_d2.A = s3.Ad;
acc_d2.B = s3.Bd;
acc_d2.C = s3.C;
acc_d2.E = s3.Ed;
acc_d2.F = s3.Kd;
acc_d2.notes = 'No bounds on process noise and observation noise are included here.';

save('data/system_examples/acc2.mat','acc_d2')