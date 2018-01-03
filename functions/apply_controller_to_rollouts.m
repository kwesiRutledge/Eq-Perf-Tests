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
%		apply_controller_to_rollouts(sys , controller , T , num_rollouts , M1 , 'missing' , missing_t )
%		apply_controller_to_rollouts(sys , controller , T , num_rollouts , M1 , 'missing' , missing_t , 'rollout_length' , rollout_length )
%
%	Inputs:
%		sys - 			This is a struct containing information about the time horizon of the model, as well as
%				  		some information on bounds of the process and measurement noises.
%						Members: .A,.B,.C,.E,.m,.d,
%
%		controller -	This struct contains the (F,u0) matrices that fully define the controller.
%						Members: .F,.u0
%
%		T - 			Time horizon for each of the rollouts
%		
%		num_rollouts - 	The number of rollouts considered.
%
%		M1 - 			Contains the important bounds for each of these controllers
%						Members: .M1
%						Optional Members: .M2
%
%		string_comm - 	A string indicating that the user would like to use a "special" option.
%						Options: 'missing','rollout_length'
%
%
%		missing_t - 	The locations of missing observations (in time) i.e. if missing_t = [1 2],
%						then the observation (or y(t)) will be not available at t = 1 and t = 2.
%
%		rollout_len -	Defines the length of each rollout. When this is not defined, the rollout length is T.
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

if (nargin < 5) or ( (nargin > 5) & (mod(nargin-5,2) == 1) )
	error(['Improper number of arguments given. (Received ' num2str(nargin) ')'])
end

sys 		= varargin{1};
controller 	= varargin{2};
T 			= varargin{3};
num_rollouts= varargin{4};
M1 			= varargin{5};

if nargin > 5
	for str_ind = 1 : mod(nargin-5,2)
		%Missing Data
		%If we assume that some points are "periodically" missing
		if strcmp(varargin{5 + 1 + (str_ind-1)*2},'missing')
			missing_t = varargin{5+str_ind*2};
		end
		%Changed Rollout Length
		%IF we assume that the rollout length IS NOT just T
		if strcmp(varargin{5+1+(str_ind-1)*2},'rollout_length')
			rollout_len = varargin{5+str_ind*2};
		end
	end
end

%%%%%%%%%%%%%%%
%% Constants %%
%%%%%%%%%%%%%%%
A_col = [];
for  i = 0:T
	A_col = [ A_col ; sys.A^i ];
end 

n = size(sys.A,1);

%Select matrix
select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

%Save important dimensions
n = size(sys.A,1);
p = size(sys.C,1);
m = size(sys.B,2);

F = controller.F;
u0 = controller.u0;

if ~isfield(sys,'E')
	sys.E = eye(n);	%If E is not explicitly defined, assume that it is identity.
end

E_bar = [];
for i = 1:T
	E_bar = blkdiag(E_bar,sys.E);
end
wd = size(sys.E,2);

% Rollout Length
if ~exist('rollout_len')
	rollout_len = T;
end

%If rollout length is defined, then use it to create these matrices.
if ~exist('missing_t')
	[S,H,Cm,~] = create_skaf_n_boyd_matrices(sys,T);
else
	[S,H,Cm,~] = create_skaf_n_boyd_matrices(sys,T,'missing',missing_t);
end

%If rollout_length is defined, then use

% disp(['size(F) = ' num2str(size(F)) ])
% disp(['size(u0) = ' num2str(size(u0)) ])
% disp(['size(H) = ' num2str(size(H))])
% disp(['size(Cm) = ' num2str(size(Cm)) ])

% Create some trajectory matrices
Pxd = (eye(n*(T+1))+H*F*inv(eye(p*T)-Cm*H*F)*Cm)*S;
Pxm = H*F*inv(eye(p*T)-Cm*H*F);
xi_factor = H*F*inv(eye(p*T)-Cm*H*F )*Cm;

% Create Disturbances
delta 	= unifrnd(-sys.d,sys.d,wd*T,num_rollouts);
mu 		= unifrnd(-sys.m,sys.m,p*T,num_rollouts);
xi_0 	= unifrnd(-M1,M1,n,num_rollouts);

if exist('missing_t')
	for missing_ob_num = missing_t
		%Remove mu for the missing times (set mu to zero)
		mu([missing_ob_num*p+1:(missing_ob_num+1)*p],:) = 0;
	end
end

%%%%%%%%%%%%%%%%%%%%%
%% Create Rollouts %%
%%%%%%%%%%%%%%%%%%%%%

xi =  A_col*xi_0 + H*repmat(u0,1,num_rollouts) + ...
        xi_factor*(A_col*xi_0 + H*repmat(u0,1,num_rollouts)) + ...
        Pxd * E_bar * delta + Pxm * mu;

xi_mag = [];
if nargout == 2
	for rollout_i = 1 : num_rollouts
		for t = 0:T
			xi_mag(t+1,rollout_i) = norm(select_m(t,T)*xi(:,rollout_i),Inf);
		end
	end
end

%%%%%%%%%%%%%%%%%%%%
%% Create Outputs %%
%%%%%%%%%%%%%%%%%%%%

varargout{1} = xi;
varargout{2} = xi_mag;

end

