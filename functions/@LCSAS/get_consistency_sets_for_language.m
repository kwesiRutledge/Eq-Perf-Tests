function [ ConsistencySets , PhiSets ] = get_consistency_sets_for_language(varargin)
%Description:
%
%Usage:
%	[ ConsistencySets , PhiSets ] = lcsas.get_consistency_sets_for_language(t,L,P_u,P_x0)
%	[ ConsistencySets , PhiSets ] = lcsas.get_consistency_sets_for_language(t,L,P_u,P_x0,'debug_flag',1)
%	[ ConsistencySets , PhiSets ] = lcsas.get_consistency_sets_for_language(t,L,P_u,P_x0,'fb_method','output')
%	[ ConsistencySets , PhiSets ] = lcsas.get_consistency_sets_for_language(t,L,P_u,P_x0,'use_proj',false)
%	[ ConsistencySets , PhiSets ] = lcsas.get_consistency_sets_for_language(t,L,P_u,P_x0,'use_proj',false,'ConsistencySetVersion',1)

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
			case 'ConsistencySetVersion'
				ConsistencySetVersion = varargin{varargin_idx+1};
				varargin_idx = varargin_idx + 2;
			otherwise
				error('Unexpected additional input.')
		end
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	if ~exist('fb_method')
		fb_method = 'output';
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

	if ~exist('ConsistencySetVersion')
		ConsistencySetVersion = 1;
	end

	n_x = P_x0.Dim;
	n_y = size(lcsas.Dyn(1).C,1);
	n_v = size(lcsas.Dyn(1).C_v,2);
	n_u = size(lcsas.Dyn(1).B,2);

	q_x0 = size(P_x0.A,1);

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	t;
	ConsistencySetVersion;

	if t > 0
		for word_idx = 1:L.cardinality()
			if debug_flag > 1
				disp(['  + word_idx = ' num2str(word_idx)]);
			end
			
			temp_lang = Language( L.words{word_idx} );

			switch ConsistencySetVersion
				case 1
					[ ConsistencySets(word_idx) , PhiSets(word_idx) ] = lcsas.consistent_set(	t,temp_lang,P_u,P_x0, ...
																								'fb_method',fb_method,'use_proj',use_proj, ...
																								'reduce_representation',reduce_flag );
				case 2
					[ ConsistencySets(word_idx) , PhiSets(word_idx) ] = lcsas.consistency_set2(	t,temp_lang,P_u,P_x0, ...
																								'fb_method',fb_method,'use_proj',use_proj, ...
																								'reduce_representation',reduce_flag );

				otherwise
					error(['Unsupported ConsistencySetVersion value: ' num2str(ConsistencySetVersion)])
			end
			

		end
	elseif t == 0
		% switch ConsistencySetVersion
		% 	case case_expression
		% 		body
		% 	otherwise
		% 		body
		% end
		for word_idx = 1:L.cardinality()
			if debug_flag > 1
				disp(['  + word_idx = ' num2str(word_idx)]);
			end
			

			first_word = L.words{word_idx};
			first_symb = first_word(1);

			q_v  = size(lcsas.Dyn(first_symb).P_v.A,1);
			temp_dyn = lcsas.Dyn(first_symb);

			if use_proj
				ConsistencySets(word_idx) = temp_dyn.C * P_x0 + temp_dyn.P_v;
			else
				ConsistencySets(word_idx) = Polyhedron();
			end

			switch ConsistencySetVersion
				case 1
					PhiSets(word_idx) = Polyhedron(	'A',[ zeros(q_x0,n_y) , P_x0.A , zeros(q_x0,n_v+n_x);
												  zeros(q_v,n_y+n_x) , temp_dyn.P_v.A , zeros(q_v,n_x);
												  zeros(q_x0,n_y+n_x+n_v), P_x0.A ], ...
											'b',[ P_x0.b ; temp_dyn.P_v.b; P_x0.b ] ,...
											'Ae',[ -eye(n_y) , temp_dyn.C , temp_dyn.C_v , zeros(n_y,n_x)  ] , ...
											'be', zeros(n_y,1) );
				case 2
					PhiSets(word_idx) = Polyhedron( 'A',[ 	zeros(q_x0,n_y) , P_x0.A ;
														temp_dyn.P_v.A*(eye(n_y)), -temp_dyn.P_v.A*temp_dyn.C ], ...
													'b', [P_x0.b; temp_dyn.P_v.b] )
				otherwise
					error(['Invalid ConsistencySetVersion: ' num2str(ConsistencySetVersion) ])
			end

		end

	else
		error(['Unexpected time input: ' num2str(t) ])

	end

end