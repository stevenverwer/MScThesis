classdef FlybyTrajectoryModel3D < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        scenario
        spacecraft
        system
        target
        T
        tStart
        tEnd
        etStart
        etEnd
    end

    methods
        function obj = FlybyTrajectoryModel3D( scenario )
            
            obj.scenario = scenario;
            obj.system = scenario.system;
            obj.target = scenario.target;

            % get spacecraft properties from scenario
            m       = scenario.satellite.mass;
            m_dry   = scenario.satellite.dryMass;
            Tmax    = scenario.propulsionSystem.thrust;
            Ve      = scenario.propulsionSystem.ve;
            
            % normalize for spacecraft & normalized properties
            obj.spacecraft.Tmax     = Tmax / m;
            obj.spacecraft.Pprop    = scenario.propulsionSystem.Pprop;
            obj.spacecraft.P0       = scenario.satellite.P0;
            obj.spacecraft.m        = m / m;
            obj.spacecraft.m_dry    = m_dry / m;
            obj.spacecraft.Ve       = Ve;
            obj.spacecraft.epsilon  = scenario.satellite.epsilon;
            obj.spacecraft.s0       = scenario.earthEscape.sf;
        end

        function setEtTime(obj, T)
            obj.tStart      = obj.target.approachDate - years(T);
            obj.tEnd        = obj.target.approachDate;
            obj.etStart     = cspice_str2et( char( obj.tStart, obj.system.dateformat ) );
            obj.etEnd       = cspice_str2et( char( obj.tEnd, obj.system.dateformat ) );
            obj.T           = T;
        end

        function [ x, vx, m, phi_x, phi_v, phi_m, xi, vi, mi, pxi, pvi, pmi, nvars] = state2vars(~, state )
            % state indices
            xi = 1:3; vi = 4:6; mi = 7; pxi = 8:10; pvi = 11:13; pmi = 14;
            nvars = 14;
            % position state
            x       = state(xi,:);
            % velocity state
            vx      = state(vi,:);
            % mass state
            m       = state(mi,:);
            % covariable position states
            phi_x   = state(pxi,:);
            % covariable velocity states
            phi_v   = state(pvi,:);
            % covariable mass states
            phi_m   = state(pmi,:);
        end

        function ds = f( obj, t, state )
            % spacecraft
            Tmax            = obj.spacecraft.Tmax;
            Ve              = obj.spacecraft.Ve;
            eps             = obj.spacecraft.epsilon;

            % time
            et = obj.etStart + (t/obj.T) * (obj.etEnd - obj.etStart);

            %x y z vx vy vz m & covariables (or lagrange multipliers)
            [ x, vx, m, phi_x, phi_v, phi_m, xi, vi, mi, pxi, pvi, pmi, nvars] = obj.state2vars( state );
            
            % ds will be the form:
            % ds(s) = A(s)*s + B(s)*u(s) + v(s,t)
            % where A is the matrix that maps the state on the derivatives
            % of the state excluding external forces; B is the matrix that
            % maps the thrust magnitude u(s) on the state and v includes
            % all the pertubations as a function of the state and time.
            
            % get bodies states and calculate influence of gravitational
            % force on state (and put in column vector v)
            dot_vx      = zeros( size( vx ) );      % influence of gravitational force on vx
            dot_phi_x   = zeros( size( phi_x ) );   % influence of gravitational force on phi_x

            for k = 2%1:length(obj.system.bodies)
                % get body properties & state
                mu      = obj.system.bodies{k}.mu;
                %[s_body, ~] = cspice_spkezr(obj.system.bodies{k}.name, et,...
                %    obj.system.frame, obj.system.abcorr, obj.system.observer);
                r_body  = [0;0;0];%conversion( s_body( 1:2 ) ,'km' ,'au' );

                % get relative position vector
                r_rel   = x - r_body;
                r       = norm( r_rel );
                
                % get accelerations caused by gravitation
                dot_vx = dot_vx - r_rel * mu / r^3;
                
                % change in the covariables
                dot_phi_x = dot_phi_x + ( mu / r^3 ) * phi_v -...
                    ( 3 * mu / r^5 ) * dot( r_rel, phi_v ) * r_rel;
            end

            % A matrix
            A               = zeros(nvars, nvars);
            A(xi, vi)       = eye(3);
            A(pvi, pxi)     = -eye(3);

            % control influence on state
            [ alpha, u ]    = obj.control_input( state );
            B               = zeros( nvars, 1 );
            B(vi)           = alpha * ( Tmax / m );
            B(mi)           = - Tmax / Ve;
            B(pmi)          = - ( Tmax / m^2 ) * norm( phi_v );

            % influence of perturbations on state derivative
            w               = zeros( nvars, 1 );
            w(vi)           = dot_vx;
            w(pxi)          = dot_phi_x;

            % calculate full state derivative
            ds              = A * state + B * u + w;
        end

        function [ alpha, u ] = control_input( obj, state )
            % get constants
            eps     = obj.spacecraft.epsilon;
            Ve      = obj.spacecraft.Ve;
            
            % get variables from state vector s
            [~,~,m,~,phi_v,phi_m] = obj.state2vars(state);
            
            % calculate optimal control direction
            alpha   = - phi_v / norm(phi_v);
            
            % calculate optimal control magnitude
            SF = norm(phi_v)* Ve / m + phi_m;
            u = 1 / ( 1 + exp( (1 - SF ) / eps ) );
        end

        function [position, isterminal, direction] = runOutOfFuelEvent(obj,t,state)
            [ ~, ~, m] = obj.state2vars( state );
            event = m - .2; % stop when running out of mass to improve solvability
            position    = [ event ]; % is positive if earth escape is possible
            isterminal  = [ 1 ];  % Halt integration
            direction   = [ 0];   % The zero can be approached from either direction
        end

        function [T, Angle, phi_x, phi_vx, phi_m] = normParams2Params(obj, input)
            T              = input(1);
            Angle          = 0;

            % for hyper unit sphere
            % coordinates to hyper unit cartesian
            span            = pi * [1 1 1 1 1 2]; 
            offset          = 0 * [ .5 .5 .5 .5 .5 .5 ]';
            beta            = diag( span ) * ( input(2:7)' - offset );
            [phi_x, phi_vx, phi_m] = obj.hypSphereAngles2UnitVec(beta,1);
        end

        function [phi_x, phi_v, phi_m] = hypSphereAngles2UnitVec(~,beta,r)
            % Basically, the langrange parameters lie on a hyper unit
            % sphere (r=1).
            % Using hyper angles to define the parameters and then mapping
            % to hyper cartesian coordinates reduces the design space by
            % one parameter and constrains the design parameters to a range
            % of 0 to 1 instead of [-inf, +inf].
            phi_x = zeros(3,1);
            phi_v = zeros(3,1);
            phi_m = zeros(1,1);
            phi_x(1) = r*cos(beta(1));
            phi_x(2) = r*sin(beta(1)) * cos(beta(2));
            phi_x(3) = r*sin(beta(1)) * sin(beta(2)) * cos(beta(3));
            phi_v(1) = r*sin(beta(1)) * sin(beta(2)) * sin(beta(3)) * cos(beta(4));
            phi_v(2) = r*sin(beta(1)) * sin(beta(2)) * sin(beta(3)) * sin(beta(4)) * cos(beta(5));
            phi_v(3) = r*sin(beta(1)) * sin(beta(2)) * sin(beta(3)) * sin(beta(4)) * sin(beta(5)) * cos(beta(6));
            phi_m(1) = r*sin(beta(1)) * sin(beta(2)) * sin(beta(3)) * sin(beta(4)) * sin(beta(5)) * sin(beta(6));
        end

        function [cost, t, y, alpha, u] = OCP( obj, input, TPBVP, forPlotting )
            % load cspice (must be loaded in order to convert time and get
            % Earth position
            cspice_furnsh( obj.system.kernels )
            
            % get converted design parameters from design parameters input
            [T, Angle, phi_x, phi_vx, phi_m] =...
                obj.normParams2Params( input );
            
            % set starting data/time
            obj.setEtTime(T);
            
            % if phi_vx==0 then 1 / norm(phi_vx) results in inf, therefore:
            if norm(phi_vx)==0
                cost = nan;
                return;
            end
            
            % index of variables in state
            % state indices
            xi = 1:3; vi = 4:6; mi = 7; pxi = 8:10; pvi = 11:13; pmi = 14;
            nvars = 14;

            % cspice get initial pos earth at etStart & target pos at etEnd
            [s_earth, ~] = cspice_spkezr(...
                obj.system.bodies{1}.name,...
                obj.etStart,...
                obj.system.frame,...
                obj.system.abcorr,...
                obj.system.observer);
            r_earth = conversion(s_earth(1:3),'km', 'au');
            v_earth = conversion(s_earth(4:6),'km/s', 'au/y');
            
            s_target = obj.target.state;
            r_target = s_target(1:3);
            v_target = s_target(4:6);
            
            % initial state and covars
            s0      = zeros(nvars,1);
            s0(xi)  = r_earth;
            s0(vi)  = v_earth;
            s0(mi)  = obj.spacecraft.s0(5);
            s0(pxi)	= phi_x;
            s0(pvi)	= phi_vx;
            s0(pmi)	= phi_m;

            % ode solver settings (some settings are for solving, plotting
            % or TPBVP using gradient-based solver
            if forPlotting % plotting high accuracy
                setOdef1 = odeset('RelTol',3e-14,'AbsTol',3e-14,...
                    'Events',@(t,s) obj.runOutOfFuelEvent(t,s));
            elseif TPBVP % gradient-based solver
                setOdef1 = odeset('RelTol',1e-13,'AbsTol',1e-13,...
                    'Events',@(t,s) obj.runOutOfFuelEvent(t,s));
            else % pso swarm accuracy
                setOdef1 = odeset('RelTol',1e-3,'AbsTol',1e-3,...
                    'Events',@(t,s) obj.runOutOfFuelEvent(t,s));
            end


            % solve the ode
            sol = ode89( @(t,s) obj.f(t,s), [0 obj.T], s0, setOdef1 );

            % for plotting
            t       = []; y       = []; alpha   = []; u       = [];
            if forPlotting
                t = [t, sol.x]; y = [y, sol.y];
                for i = 1:length(sol.x)
                    [alphat, ut] = obj.control_input(sol.y(:,i));
                    alpha = [alpha, alphat]; u     = [u, ut];
                end
            end
            
            % check end conditions (if the spacecraft ran out of
            % propellant)
            if isempty(sol.ie)
            elseif any(sol.ie==1)
                cost = nan;
                cspice_kclear;
                return
            elseif any(sol.ie==2)
                cost = nan;
                cspice_kclear;
                return
            end
            
            % final state
            s_tf = sol.y(:,end);

            % boundary constraints
             ceq = [
                s_tf(xi) - r_target;
                .1 * s_tf(pvi); % weighted lagrange parameters constraints to improve solvability
                .1 * s_tf(pmi); % weighted lagrange parameters constraints to improve solvability
                ];
            if TPBVP
                cost = ceq;
            else
                cost = sum(ceq.^2); % squared sum of boundary contraints
            end
            
            % close cspice
            cspice_kclear;
        end

        function stop = pswplotranges(obj,optimValues,state)
            stop = false; % This function does not stop the solver
            switch state
                case 'init'

                case 'iter'
                    xlim([-1.1 1.1])
                    ylim([-1.1 1.1])
                    grid minor
                    xlabel('$AU$','Interpreter','latex')
                    ylabel('$AU$','Interpreter','latex')
                    title('Cruise from Earth to asteroid current best guess','Interpreter','latex')
                    [fval, t, y, alpha, u] = obj.OCP( optimValues.bestx, true, true);
                    yb = zeros(size(y));
                    yb(1:3,:) = conversion(y(1:3,:), 'au', 'km');
                    yb(4:6,:) = conversion(y(4:6,:), 'au/y', 'km/s');
                    yb(7:14,:) = y(7:14,:);
                    cspice_furnsh( obj.system.kernels )
                    time = t;
                    obj.system.sol_ets                  = zeros(1,length(time));
                    for n=1:length(time)
                        obj.system.sol_ets(n)           = cspice_str2et(char(obj.tStart + years(time(n))));
                    end
                    
                    [earth.sol_states, ~] = cspice_spkezr(...
                        obj.system.bodies{1}.name,...
                        obj.system.sol_ets,...
                        obj.system.frame,...
                        obj.system.abcorr,...
                        obj.system.observer);
                    
                    yb_earth = yb(1:3,:) - earth.sol_states(1:3,:);
                    obj.target.sol_states_earth = obj.target.state - earth.sol_states;
                    drawnow % Show the plot
                    obj.target.sol_states_earth = obj.target.state - earth.sol_states;

                    plot(0,0,'o','Color',[0.9290 0.6940 0.1250],'MarkerFaceColor',[0.9290 0.6940 0.1250],'MarkerSize',10)
                    hold on
                    plot(  conversion(earth.sol_states(1,1),'km','au'),...
                            conversion(earth.sol_states(2,1),'km','au'),'x','Color',[0 0.4470 0.7410],'LineWidth',1)
                    plot(  conversion(earth.sol_states(1,:),'km','au'),...
                            conversion(earth.sol_states(2,:),'km','au'),'--','Color',[0 0.4470 0.7410],'LineWidth',1)
                    plot(  conversion(earth.sol_states(1,end),'km','au'),...
                            conversion(earth.sol_states(2,end),'km','au'),'o','Color',[0 0.4470 0.7410],'MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410])
                    plot(obj.target.state(1,:),...
                          obj.target.state(2,:),'v','Color','black','LineWidth',1)%,'MarkerFaceColor','black')
                    
                    plot(y(1,1),y(2,1),'.','Color','red')
                    plot(y(1,:),y(2,:),'-','Color','red','LineWidth',1)
                    hold off
                    cspice_kclear;
                case 'done'
                    % No cleanup necessary
            end
        end
    end

end