function [ ConsistencySets , PhiSets ] = get_consistency_sets_for_language(varargin)
%Description:
%
%Usage:
%	[ ConsistencySets , PhiSets ] = lcsas.get_consistency_sets_for_language(t,L,P_u,P_x0)
%	[ ConsistencySets , PhiSets ] = lcsas.get_consistency_sets_for_language(t,L,P_u,P_x0,'debug_flag',1)
%	[ ConsistencySets , PhiSets ] = lcsas.get_consistency_sets_for_language(t,L,P_u,P_x0,'fb_method','output')

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if nargin < 5
		error('Not enough input arguments.')
	end

	lcsas = varargin{1};
	t = varargin{2};
	L = varargin{3};
	P_u = varargin{4};
	P_x0 = varargin{5};

	if ~isa(lcsas,'LCSAS')
		error('Expecting the first input to be a LCSAS object.')
	end

	if ~isa(L,'Language')
		error('Expecting the language input to be a Language object.')
	end

	varargin_idx = 6;
	while varargin_idx <= nargin
		switch varargin{varargin_idx}
			case 'fb_method'
				fb_method = varargin{varargin_idx+1};
				if ~(strcmp(fb_method,'state') || strcmp(fb_method,'output'))
					error(['Invalid feedback type: ' fb_method ])
				end
				varargin_idx = varargin_idx + 2;
			case 'use_proj'
				use_proj = varargin{varargin_idx+1};
				if ~islogical( use_proj )
					error('The flag for ''use_proj'' should be a boolean.')
				end
				varargin_idx = varargin_idx + 2;
			case 'reduce_representation'
				reduce_flag = varargin{varargin_idx+1};
				if ~islogical( reduce_flag )
					error('The flag for ''reduce_flag'' should be a boolean.')
				end
				varargin_idx = varargin_idx + 2;
			case 'debug_flag'
				debug_flag = varargin{varargin_idx+1};
				varargin_idx = varargin_idx + 2;
			otherwise
				error('Unexpected additional input.')
		end
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	if ~exist('fb_method')
		fb_method = 'state';
	end

	if ~exist('use_proj')
		use_proj = true;
	end

	if ~exist('reduce_flag')
		reduce_flag = true;
	end

	if ~exist('debug_flag')
		debug_flag = 0;
	end

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	% PhiSets = [];
	% ConsistencySets = [];

	for word_idx = 1:L.cardinality()
		if debug_flag > 1
			disp(['  + word_idx = ' num2str(word_idx)]);
		end
		
		temp_lang = Language( L.words{word_idx} );

		[ ConsistencySets(word_idx) , PhiSets(word_idx) ] = lcsas.consistent_set(t,temp_lang,P_u,P_x0,'fb_method',fb_method);

	end

end