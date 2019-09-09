function [ fb_law , opt_out ] = rec_synthesis(varargin)
%Description:
%	Synthesizes a controller that verifies if the system can guarantee a return to an initial set X0 after the switching behavior
%	of the input language is done.
%
%Usage:
%	[fb_law,opt_out] = lcsas.rec_synthesis(L,X0)
%	[fb_law,opt_out] = lcsas.rec_synthesis(L,X0,'Pu',Pu)


	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if nargin < 3
		error('Expected at least 3 inputs. Please see usage comment. (help LCSAS.rec_synthesis)')
	end

	lcsas = varargin{1};
	L = varargin{2};
	X0 = varargin{3};

	argin_idx = 4;
	while argin_idx <= nargin
		switch varargin{argin_idx}
			case 'Pu'
				P_u = varargin{argin_idx+1};
				argin_idx = argin_idx + 2;
			case 'verbosity'
				verbosity = varargin{argin_idx};
				argin_idx = argin_idx + 2;
			otherwise
				error(['Unexpected input to the function!: ' varargin{argin_idx} ])
		end
	end

	if ~exist('verbosity')
		verbosity = 1;
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	n = size(lcsas.Dyn(1).A,1);
	m = size(lcsas.Dyn(1).B,2);
	p = size(lcsas.Dyn(1).C,1);
	wd = size(lcsas.Dyn(1).B_w,2);
	vd = size(lcsas.Dyn(1).C_v,2);

	select_m = @(t,T_r) [zeros(n,t*n), eye(n), zeros(n,(T_r-t)*n) ];

	ops = sdpsettings('verbose',verbosity);

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	disp('Hi!!')

end