%lead_foll_tests.m


%% Test 1: Verifying that offset calculation works

cube_dimx = 2;
cube_dimy = 2;
dt = 0.1;

[cont_ad,disc_ad,offset_vals] = get_lead_follow_aff_dyn2(cube_dimx,cube_dimy,dt);

assert(isequal(offset_vals(:,1),[-2;1]));

%% Test 2: Verifying that the origin is really an equilibria for the dynamical system

cube_dimx = 2;
cube_dimy = 2;
dt = 0.1;

foll_cube.n = cube_dimx*cube_dimy;

[cont_ad,disc_ad,offset_vals] = get_lead_follow_aff_dyn2(cube_dimx,cube_dimy,dt);

%Attempt to plot trajectories for a given initial condition
u0 = zeros(2,1);
ode = @(t,x) cont_ad.A*x+cont_ad.B*u0 + cont_ad.f;

[t,y] = ode45(ode,[0,10],ones(1,2*(foll_cube.n+1)));

for drone_ind = 1:1+foll_cube.n
	figure;
	plot(y(:,drone_ind),y(:,drone_ind+1+foll_cube.n))
end

figure;
hold on;
for drone_ind = 1:1+foll_cube.n
	plot(y(:,drone_ind),y(:,drone_ind+1+foll_cube.n))
end

figure;
subplot(2,1,1)
hold on;
for drone_ind = 1:1+foll_cube.n
	plot(t,y(:,drone_ind))
end

subplot(2,1,2)
hold on;
for drone_ind = 1:1+foll_cube.n
	plot(t,y(:,drone_ind+1+foll_cube.n))
end
