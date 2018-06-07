function [  ] = filter_synthesis( varargin )
	%Description:
	%	Creates a filter which should obey the specification that:
	%	- For all times at which the state of the automaton s(t) is a state that contains the recovery node,
	%	  the estimation error ||\xi(t)|| <= M1
	%	- For all other times t, ||\xi(t)|| <= M2.
	%
	%Examples:
	%	[  ] = oa.filter_synthesis( ad , M1 , M2 )
	%	[  ] = oa.filter_synthesis( ad , M1 , M2 , 'verbosity' , verbosity )

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if nargin < 4
		error(['There are not enough input arguments! Only ' num2str(nargin) ' given.'])
	end

	oa = varargin{1};
	ad = varargin{2};
	M1 = varargin{3};
	M2 = varargin{4};

	verbosity = 2;

	if nargin > 4
		str_in = varargin{5};
		switch str_in
		case 'verbosity'
			verbosity = 1;
		end
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	n = size(ad.A,1);
	m = size(ad.B,2);
	p = size(ad.C,1);
	wd = size(ad.B_w,2);
	vd = size(ad.C_v,2);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Organize the Paths that We're Interested in %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	[L,S] = oa.find_all_2R_paths();

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Perform Synthesis for Every Subset of Words in L that Start with the Same State S %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	S_covered = [];

	S_queue = [ oa.OX_0 ];
	
	%
	L0 = {}; S0 = {};
	for S_ind = 1:length(S)
		if S{S_ind}(1) == oa.OX_0 
			L0{S_ind} = L{S_ind};
			S0{S_ind} = S{S_ind};
		end
	end
	%Update the interesting start state queue.
	for S_ind = 1:length(S0)
		%
		if ~any(S0{S_ind}(end) == S_covered )
			%If the end of the state sequence is not already considered, then
			%add it to the states tp cpmsodered om S_queue.
			S_queue = [ S_queue ; S0{S_ind}(end) ];
		end
	end
	S_queue = S_queue(2:end);
	
	%Apply Synthesis
	[ opt_data , contr ] = free_rec_design_pb( ad , 'Min_M3' , M1 , M2 , L0 );

	while( ~isempty(S_queue) )
		%Select State Index from Queue
		s0 = S_queue(1);

		

	end

	disp('Finished filter_synthesis.')

end