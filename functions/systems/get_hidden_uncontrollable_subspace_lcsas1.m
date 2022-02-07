function [ lcsas_out , P_u , Pw1 , x0 , P_target ] = get_hidden_uncontrollable_subspace_lcsas1(varargin)
        %Description:
        %       Constructs Liren Yang's uncontrollable subspace example from his Zonotope paper.
        %       The subspace has 4d in his example, but this one should have a defined size.
        %Usage:
        %       [ lcsas_out , P_u , Pw1 , x0 , P_target ] = get_hidden_uncontrollable_subspace_lcsas1('UncontrollableSubspaceDimension',4)

        % x+ = Ax + Bu + Ew + K
        %  y = Cx + Fv

        % Input Processing

        [ lcsas_settings , TimeHorizon , dt , us_dim ] = get_hidden_uncontrollable_subspace_lcsas1_ip( varargin{:} );

        % Constants
        n_x = 3*2 + us_dim;
        n_y = n_x;
        n_u = 3;

        % Create System Matrices

        A_ct = zeros(n_x);
        A_ct([1:3*2],[1:3*2]) = [0,1,0,0,0,0;
                                 0,0,0,0,0,0;
                                 0,0,0,1,0,0;
                                 0,0,0,0,0,0;
                                 0,0,0,0,0,1;
                                 0,0,0,0,0,0];

        for repeat_index = 1:(us_dim/2)
                A_ct(6+(repeat_index-1)*2+[1:2],6+(repeat_index-1)*2+[1:2]) = [ -10^(-2*repeat_index) , repeat_index ; -repeat_index , -10^(-2*repeat_index) ];
        end

        for j_index = 7:2:6+us_dim
                
                switch mod(j_index-1,6)
                case 0
                        %First entry
                        A_ct([1:n_x-us_dim],j_index+[0:1]) = [  1 , 0 ;
                                                                0 , 0 ;
                                                                0 , -1;
                                                                0 , 0 ;
                                                                0 , 0 ;
                                                                0 , 0 ];
                case 2
                        %Second part of pattern
                        A_ct([1:n_x-us_dim],j_index+[0:1]) = [  0 , 1 ;
                                                                0 , 0 ;
                                                                0 , 0 ;
                                                                0 , 0 ;
                                                                1 , 0 ;
                                                                0 , 0 ];

                case 4
                        %Third part of pattern
                        A_ct([1:n_x-us_dim],j_index+[0:1]) = [  0 , 0 ;
                                                                0 , 0 ;
                                                                -1 , 0;
                                                                0 , 0 ;
                                                                0 , 1 ;
                                                                0 , 0 ];
                otherwise
                        error(['There was an issue writing down the A matrix!'])
                end

        end

        % A_ct = [0,1,0,0,0,0,1,   0,    0,     1;
        %         0,0,0,0,0,0,0,      0,    0,     0;
        %         0,0,0,1,0,0,0,     -1, 0,     0;
        %         0,0,0,0,0,0,0,      0,    0,     0;
        %         0,0,0,0,0,1,0,      0,    1,   0;
        %         0,0,0,0,0,0,0,      0,    0,     0;
        %         0,0,0,0,0,0,-1e-2,  1,    0,     0;
        %         0,0,0,0,0,0,-1,    -1e-2, 0,     0;
        %         0,0,0,0,0,0,0,      0,   -1e-4,  2;
        %         0,0,0,0,0,0,0,      0,   -2,     -1e-4];


        B_ct = zeros(n_x,n_u);
        B_ct([1:6],:) = [0 0 0;
                        1 0 0;
                        0 0 0;
                        0 1 0;
                        0 0 0;
                        0 0 1];

        % B_ct = [0 0 0;
        %         1 0 0;
        %         0 0 0;
        %         0 1 0;
        %         0 0 0;
        %         0 0 1;
        %         0 0 0;
        %         0 0 0;
        %         0 0 0;
        %         0 0 0];

        K_ct = zeros(n_x,1);
        E_ct = eye(n_x);

        sys_ct = ss(A_ct,B_ct,eye(n_x),0);
        sys_dt = c2d(sys_ct, dt,'zoh');
        A = sys_dt.A;
        B = sys_dt.B;

        sys_ct = ss(A_ct,E_ct,eye(n_x),0);
        sys_dt = c2d(sys_ct, dt,'zoh');
        E = sys_dt.B;

        sys_ct = ss(A_ct,K_ct,eye(n_x),0);
        sys_dt = c2d(sys_ct, dt,'zoh');
        K = sys_dt.B;

        scl = lcsas_settings.scl;
        Pw1 = Polyhedron('lb',scl*[-0.12,-0.2, -0.12,-0.2, -0.08,-0.2, -0.1*ones(1,us_dim)],'ub',scl*[0.12, 0.2, 0.12, 0.2, 0.12, 0.2,0.1*ones(1,us_dim)]);

        ad1 = Aff_Dyn( ...
                A,B,K,zeros(n_y), ...
                Polyhedron('lb',-1*ones(1,n_y),'ub',1*ones(1,n_y)),Pw1 );

        ad2 = ad1;
        ad2.A(4,2) = 0.2;


        x0 = [ -60 ; 0 ; -60 ; 0 ; -70 ; 0 ; scl*ones(us_dim,1) ];
        P_x0 = Polyhedron('lb',x0','ub',x0');

        %% Create Outputs

        P_u = Polyhedron('lb',-lcsas_settings.eta_u*ones(1,n_u),'ub',lcsas_settings.eta_u*ones(1,n_u));

        lcsas_out = LCSAS( ...
                [ad1,ad2] , Language(1*ones(1,TimeHorizon),2*ones(1,TimeHorizon)) , ...
                'X0' , P_x0, ...
                'U' , P_u );

        P_target = lcsas_settings.P_target;

end

function [ lcsas_settings , TimeHorizon , dt , us_dim ] = get_hidden_uncontrollable_subspace_lcsas1_ip( varargin )
        %Description:
        %       Parses the potential nonstandard inputs to the lcsas.

        %% Default

        lcsas_settings = struct( ...
                'scl', 0.1 , ...
                'dt', 0.5, ...
                'UncontrollableSubspaceDimension', 4, ...
                'eta_u', 10 , ...
                'TimeHorizon', 5 , ...
                'P_target', []  ...
                );

        %% Input Processing %%

        argument_index = 1;
        while nargin >= argument_index
                switch varargin{argument_index}
                case 'TimeHorizon'
                        lcsas_settings.TimeHorizon = varargin{argument_index+1};
                case 'dt'
                        lcsas_settings.dt = varargin{argument_index+1};
                case 'scl'
                        lcsas_settings.scl = varargin{argument_index+1};
                case 'UncontrollableSubspaceDimension'
                        lcsas_settings.UncontrollableSubspaceDimension = varargin{argument_index+1};
                case 'P_target'
                        lcsas_settings.P_target = varargin{argument_index+1};
                case 'eta_u'
                        lcsas_settings.eta_u = varargin{argument_index+1};
                otherwise
                        error(['Unexpeced input to get_uncertain_thick_pendulum_lcsas: ' varargin{argument_index} ])
                end
                % Increment
                argument_index = argument_index + 2;
                
        end

        %% Checking Inputs
        if mod(lcsas_settings.UncontrollableSubspaceDimension,2) ~= 0
                error(['The uncontrollable subspace dimension must be a multiple of 2, received ' num2str(lcsas.UncontrollableSubspaceDimension) '.'] )
        end

        %% Creating Some Outputs %%

        TimeHorizon = lcsas_settings.TimeHorizon;
        dt          = lcsas_settings.dt;
        us_dim      = lcsas_settings.UncontrollableSubspaceDimension;

        % If P_target hasn't been created yet, then create it!
        if isempty(lcsas_settings.P_target)
                target_lb = [ -80 , -10 , -80 , -10 , -60 , -10 , -ones(1,us_dim)*10^(3) ];
                target_ub = [ -40 , 10 , -40 , 10 , -20 , 10 , ones(1,us_dim)*10^(3) ];

                % Create P_target
                lcsas_settings.P_target = Polyhedron('lb',target_lb,'ub',target_ub);

        end

end
