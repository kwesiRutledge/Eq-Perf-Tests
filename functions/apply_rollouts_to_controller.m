function [ varargout ] = apply_controller_to_rollouts(varargin)
% apply_rollouts_to_controller.m
%	Description:
%		The objective of this function is to output the result of applying the controller defined by
%		"controller" to num_rollouts of rollouts for the system "sys"
%
%	Model:
%		The dynamics for this system are assumed to be the following ()
%
%     		xi[k+1] = A xi[k] + B u[k] + w[k]
%     		y[k]    = C xi[k] + v[k]
%
%		The controller is entirely parameterized by 2 matrices (F,u0).
%		
%			u[k] = \sum_{t=0}^k F(k,t) * y[t] + u0[k]
%
%	Usage:
%		apply_controller_to_rollouts(sys , controller , T , num_rollouts , M1)
%		apply_controller_to_rollouts(sys , controller , T , num_rollouts , M1 ,'missing' , missing_t )
%
%	Inputs:
%		sys - 			This is a struct containing information about the time horizon of the model, as well as
%				  		some information on bounds of the process and measurement noises.
%						Members: .A,.B,.C,.m,.d,
%
%		controller -	This struct contains the (F,u0) matrices that fully define the controller.
%						Members: .F,.u0
%
%		T - 			Time horizon for each of the rollouts
%		
%		num_rollouts - 	The number of rollouts considered.
%
%		bound_list - 	Contains the important bounds for each of these controllers
%						Members: .M1
%						Optional Members: .M2
%
%	Outputs:
%		xi - 		Matrix containing every one of the 'num_rollouts' of rollouts
%					(each column is a rollout)
%
%		xi_mag -	Infinity_norm of each state in each rollout sequence instead
%					of the full state value	
%

%%%%%%%%%%%%%%%%%%%%%%
%% Input Processing %%
%%%%%%%%%%%%%%%%%%%%%%

if (nargin < 5)
	error(['Improper number of arguments given. (Received ' num2str(nargin) ')'])
end

sys 		= varargin{1};
controller 	= varargin{2};
T 			= varargin{3};
num_rollouts= varargin{4};
M1 			= varargin{5};

if nargin > 5
	if strcmp(varargin{6},'missing')
		missing_t = varargin{7}
	end
end

%%%%%%%%%%%%%%%
%% Constants %%
%%%%%%%%%%%%%%%
if ~exist('missing_t')
	[S,H,Cm,~] = create_skaf_n_boyd_matrices(sys,T);
else
	[S,H,Cm,~] = create_skaf_n_boyd_matrices(sys,T,'missing',missing_t);
end

A_col = [];
for  i = 0:T
	A_col = [ A_col ; sys.A^i ];
end 

%Save important dimensions
n = size(sys.A,1);
p = size(sys.C,1);
d_u = size(sys.B,2);

F = controller.F;
u0 = controller.u0;

% disp(['size(F) = ' num2str(size(F)) ])
% disp(['size(u0) = ' num2str(size(u0)) ])
% disp(['size(H) = ' num2str(size(H))])
% disp(['size(Cm) = ' num2str(size(Cm)) ])

% Create some trajectory
Pxd = (eye(n*(T+1))+H*F*inv(eye(p*T)-Cm*H*F)*Cm)*S ;
Pxm = H*F*inv(eye(p*T)-Cm*H*F);
xi_factor = (eye(n*(T+1)) +  H*F*inv(eye(p*T)-Cm*H*F )*Cm);

% 

% Create Disturbances
delta 	= unifrnd(-sys.d,sys.d,n*T,num_rollouts);
mu 		= unifrnd(-sys.m,sys.m,p*T,num_rollouts);
xi_0 	= unifrnd(-M1,M1,n,num_rollouts);

%%%%%%%%%%%%%%%%%%%%%
%% Create Rollouts %%
%%%%%%%%%%%%%%%%%%%%%

xi =   xi_factor*(A_col*xi_0 + H*u0) + Pxd * delta + Pxm * mu;

xi_mag = [];
if nargout == 2
	for rollout_i = 1 : num_rollouts
		for t = 0:T
			xi_mag(t+1,rollout_i) = norm(xi(n*t+1:n*(t+1),rollout_i),Inf);
		end
	end
end

%%%%%%%%%%%%%%%%%%%%
%% Create Outputs %%
%%%%%%%%%%%%%%%%%%%%

varargout{1} = xi;
varargout{2} = xi_mag;

end

