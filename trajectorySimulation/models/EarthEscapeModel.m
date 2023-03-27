classdef EarthEscapeModel < handle
    %EarthEscapeModel simulates the trajectory from Earth orbit until Earth
    %escape

    properties
        scenario % details about the initial orbit and other circumstances
        spacecraft % details about the spacecraft
        system % the coordinate system & dynamics as part of the scenario
    end

    methods
        function obj = EarthEscapeModel( scenario )
            obj.scenario    = scenario;
            obj.system      = scenario.system;

            % get spacecraft properties from scenario
            m               = scenario.satellite.mass;
            m_dry           = scenario.satellite.dryMass;
            Tmax            = scenario.propulsionSystem.thrust;
            Ve              = scenario.propulsionSystem.ve;
            
            % normalize for spacecraft based on spacecraft mass
            obj.spacecraft.Tmax     = Tmax / m;
            obj.spacecraft.Pprop    = scenario.propulsionSystem.Pprop;
            obj.spacecraft.P0       = scenario.satellite.P0;
            obj.spacecraft.m        = m / m;
            obj.spacecraft.m_dry    = m_dry / m;
            obj.spacecraft.Ve       = Ve;
            obj.spacecraft.Cd       = scenario.satellite.dragCoefficient;
            obj.spacecraft.A        = scenario.satellite.frontalArea;
        
            % get starting orbit properties
            obj.spacecraft.apoapsis     = scenario.satellite.apoapsis;
            obj.spacecraft.periapsis    = scenario.satellite.periapsis;
        end

        function [s0] = startingOrbit(obj)
            % get starting orbit properties
            a       = (obj.spacecraft.apoapsis + obj.spacecraft.periapsis) / 2;

            % Earth standard gravitational paramater
            mu      = obj.system.bodies{1}.mu;

            % starting distance from Earth
            r       = obj.spacecraft.periapsis;

            % starting velocity at apoapsis
            vth     = sqrt( mu * ( 2 / r - 1 / a ) );

            % starting state
            s0      = [r; 0; 0; vth; obj.spacecraft.m];
        end

        function [ r, th, vr, vth, m, ri, thi, vri, vthi, mi ] = state2vars(~, state )
            % state indices
            ri = 1; thi = 2; vri = 3; vthi = 4; mi = 5;
            % distance states
            r       = state(ri,:);
            % orbit angle states
            th      = state(thi,:);
            % velocity parallel to distance vector states
            vr       = state(vri,:);
            % velocity perpendicular to distance vector states
            vth     = state(vthi,:);
            % mass states
            m       = state(mi,:);
        end

        function [ events, isterminal, direction ] = odeEvents( obj, t, state )
            % get state variables
            [ r, th, vr, vth, m ] = obj.state2vars( state );
            
            % Earth escape event
            event1 = ( vr^2 + vth^2 )/2 - obj.system.bodies{1}.mu/r;

            % Run out of fuel event
            event2 = obj.spacecraft.m_dry - m;
            
            % Events
            events      = [ event1, event2 ];
            isterminal  = [ 1,  1 ];  % Halt integration
            direction   = [ 0,  0 ];  % The zero can be approached from either direction
        end

        function [ ds, u, alpha, drag_r, drag_th ] = f( obj, t, state )
            % get system variables
            %Earth mu
            mu          = obj.system.bodies{1}.mu;
            
            % spacecraft variables
            Tmax        = obj.spacecraft.Tmax; %1AU
            Ve          = obj.spacecraft.Ve;
            Cd          = obj.spacecraft.Cd;
            A           = obj.spacecraft.A;
        
            % get state variables
            [ r, th, vr, vth, m, ri, thi, vri, vthi, mi ] = obj.state2vars( state );
            
            % heuristic thrusting law
            alpha       = atan2( vr, vth );
            u           = 1;
            
            % ATMOSPHERIC DRAG & normalize
            h           = conversion( r - obj.system.bodies{1}.r, 'au', 'm' );
            if h>1e6
                rho             = conversion( 9e-27, '1/m^3', '1/au^3');
            else
                atmosDensity    = @(h) 10.^ppval(obj.system.atmosphericDensityMakimaPp, log10(h));
                rho             = conversion( atmosDensity(h), '1/m^3', '1/au^3' );
            end
            drag_r      = -.5 * Cd * A * rho * vr*norm(vr) * (1/obj.scenario.satellite.mass);
            drag_th     = -.5 * Cd * A * rho * vth*norm(vth) * (1/obj.scenario.satellite.mass);

            % A matrix
            A           = zeros( length(state), length(state) ); 
            A(ri,vri)   = 1;
            A(thi,vthi) = 1 / r;
            A(vri,vthi) = vth / r;
            A(vthi,vri) = - vth / r;
            
            % B matrix
            B           = zeros( size( state ) );
            B(vri)      = sin(alpha) * (Tmax / m);
            B(vthi)     = cos(alpha) * (Tmax / m);
            B(mi)       = - (Tmax / Ve);
            
            % perturbations
            w           = zeros( size( state ) );
            
            % gravity
            w(vri)      = w(vri) - mu / norm( r )^2;

            % drag
            w(vri)      = w(vri)  + drag_r / m;
            w(vthi)     = w(vthi) + drag_th / m;
            
            % state derivative
            %ds/dt = A.s + B.u + w
            ds          = A * state + B * u + w;
        end
    end
end