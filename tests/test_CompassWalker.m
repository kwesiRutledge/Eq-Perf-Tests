%test_CompassWalker.m
%Description:
%	Describes how to use the CompassWalker class


function tests = test_CompassWalker
	%disp(localfunctions)
	tests = functiontests(localfunctions);

function include_relevant_libraries()
	%Description:
	%	Attempts to add the relevant libraries/toolboxes to the path.

	%% Constants

	is_on_personal_mac = strcmp(getenv('USER'),'kwesirutledge');
	is_on_great_lakes = (strcmp(getenv('USER'),'krutledg') && isunix) ;

	addpath('../')

	%% Algorithm
	if is_on_personal_mac
		include_fcns2({'mosek','gurobi','tbxmanager'},'PathToDirectoryWithToolboxes','../../../')
	elseif is_on_great_lakes
		include_fcns2({'tbxmanager','YALMIP'},'PathToDirectoryWithToolboxes','../../../')
	else
		warning('Warning: There is no guaranteed performance for systems that are not on Kwesi''s computer.')
	end

	addpath(genpath('../functions'));

function test1_constructor(testCase)
	%Description:
	%	Tests constructor

	%% Including Libraries
	include_relevant_libraries();

	%% Constants
	cw = CompassWalker();

	%% Algorithm

	assert( (cw.a == 0.5) && (cw.b == 0.5) && (cw.m == 5) && ...
			(cw.m_H == 10) && (cw.l == 1) && (cw.phi == pi/6) )

function test1_alpha(testCase)
	%Description:
	%	Tests the subfunction alpha

	%% Including Libraries
	include_relevant_libraries();

	%% Constants
	cw = CompassWalker();
	eps1 = 10^(-3);

	%% Algorithm

	cw.CurrentState = [ 3*pi/2 ; 2*pi ; 0 ; 0 ];

	assert( abs( cw.alpha() - (pi/4) ) < eps1 )

function test1_visualize_trajectory(testCase)
	%Description:
	%	Makes sure that the error-handling of visualize_trajectory() is working properly.
	%	Submits a trajectory that does not have the proper dimensions.

	%% Including Libraries
	include_relevant_libraries();

	%% Constants
	cw = CompassWalker();
	
	temp_traj = randn(3,20);

	%% Algorithm

	try
		cw.visualize_trajectory( temp_traj )
		assert(false)
	catch e
		assert(strcmp(e.message,['The input trajectory x_trajectory should have 4 rows, instead it has 3.' ]))
	end

function test1_M(testCase)
	%Description:
	%	Computes the matrix M for the compass walker which is a function of the
	%	two walker foot angles.

	%% Including Libraries
	include_relevant_libraries();

	%% Constants
	cw1 = CompassWalker();
	cw1.CurrentState = [ 3*pi/2 ; 2*pi ; 0 ; 0 ];

	m = cw1.m;
	m_H = cw1.m_H;
	a = cw1.a;
	b = cw1.b;
	l = cw1.l;

	theta_ns = cw1.CurrentState(1);
	theta_s = cw1.CurrentState(2);

	%% Algorithm
	M1 = [ m*(b^2) , -m*l*b*cos( theta_s - theta_ns ) ;
			-m*l*b*cos(theta_s - theta_ns) , (m_H+m)*l^2+m*a^2 ];

	assert( all(all( M1 == cw1.M( cw1.CurrentState ) )) )

function test1_G(testCase)
	%Description:
	%	Computes the matrix G for the compass walker which is a function of the
	%	two walker foot angles.

	%% Including Libraries
	include_relevant_libraries();

	%% Constants
	cw1 = CompassWalker();
	cw1.CurrentState = [ 3*pi/2 ; 2*pi ; 0 ; 0 ];

	m = cw1.m;
	m_H = cw1.m_H;
	a = cw1.a;
	b = cw1.b;
	l = cw1.l;

	theta_ns = cw1.CurrentState(1);
	theta_s = cw1.CurrentState(2);

	g = 10;

	%% Algorithm

	G = [ m*b*sin(theta_ns) ; -(m_H*l+m*a+m*l)*sin(theta_s) ] * g;

	assert( all(all( G == cw1.G( cw1.CurrentState ) )) )

function test1_toLCSAS(testCase)
	%Description:
	%	Tests the generation of a 