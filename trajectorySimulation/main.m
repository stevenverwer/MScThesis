%{
This script optimizes a low thrust trajectory using 
Global Optimization Toolbox's Partical Swarm Optimization
while using an optimal control strategy using the calculus of variations
indirect single shooting method

The first part consists of an optimal earth escape trajectory based on heuristic laws
and the second part is a deep space cruise to the target location.

dependencies:
    Aerospace Toolbox
    Global Optimization Toolbox
    Parallel Processing Toolbox
    MICE
    JPL Horizons
    SPICE KERNELS
    propulsionSystems.xlsx
    satellites.xlsx
    targets.xlsx
    conversion.m

Parts:
    Cases: defines the different scenarios and combinations that should be
    simulated.

    Constants: defines the constants of the trajectory such as bodies,
    target, spacecraft and the trajectories of bodies and target.
    
    Then the initial orbit is defined as 90000 - 295 km altitude SSGTO
    around earth in the EarthEscapeModel.m class.

    Finally, the states of the satellite are defined as pos/vel/mass and
    their corresponding lagrange multipliers as used in the Hamiltonian.

    FlybyTrajectoryModel3D.m class handles the definition of all the 
    dynamics, control and optimization cost function.

    To change scenario, propulsion system, or satellite change:
     -->propulsionSystems.xlsx
     -->satellites.xlsx
     -->targets.xlsx
    For more cases, just add more rows the code will calculate every
    combination.
%}

% add paths to functions, models, data and results folder
addpath('functions');
addpath('models');
addpath('data');
addpath('results');

% NASA's SPICE directories DO NOT DELETE!
addpath('naif\mice\lib');
addpath('naif\mice\src\mice');
addpath('naif\kernels');

% settings (for NASA's SPICE) ---------------------------------------------
system.dateformat   = 'MMM dd HH:mm:ss.SSSSSSSSS yyyy';
system.frame        = 'ECLIPJ2000'; %   ECLIPJ2000 reference frame, targets 
                                % should also be in this frame.
system.abcorr       = 'None';
system.observer     = '0';
system.kernels      = 'naif\kernels\targets.furnsh'; % not necesarry for 
                                % targets, but it contains other bodies

% TIME --------------------------------------------------------------------
max_mission_duration = 5; %years

% BODIES ------------------------------------------------------------------
% body constants are:
%{
    name    = the body's name
    r       = radius of body in astronomical units
    mu      = standard gravitational parameter of a body in au^3 per year^2
%}

% Earth
earth.name          = 'Earth';
earth.r             = conversion( 6378.1363e3, 'm', 'au');
earth.mu            = conversion( 3.9860044188e14, 'm^3/s^2', 'au^3/y^2' );

% Sun
sun.name            = 'Sun';
sun.mu              = conversion( 1.327124400189e20, 'm^3/s^2', 'au^3/y^2' );
sun.r               = conversion( 696340000, 'm', 'au' );

% Moon
moon.name           = 'Moon';
moon.r              = conversion( 1737.5e3, 'm', 'au' );
moon.mu             = conversion( 4.90486959e12, 'm^3/s^2', 'au^3/y^2' );

system.bodies = {earth};
% SELECT BODIES & TARGET   ------------------------------------------------
system.bodies = [system.bodies, ...
    { sun, ...
      moon }];

% keep the workspace clean, delete what you don't need.
clear sun moon earth;

% AIRDENSITIES FOR DRAG   -------------------------------------------------
% plot airdensities:
plotAirdensities = false; % set to true to plot
%{
    to obtain the drag in Earth orbit, the NRLMSISE-00 is used. This
    particular model only works from 0 to 10^6 meter. This datapoints
    obtained from atmosnrlmsise00 is fitted with a makima interpolation
    method which produces piecewise polynomials with continuous first-order
    derivatives. The algorithm avoids excessive local undulations.

    variables:
        h       --> height
        rhos    --> densities of atmospheric components
        rho     --> air density
%}

% F10.7 average (default = 150)
F107a = 300;

% F10.7 daily   (default = 150)
F107 = 300;

% magnetic index (default = 4)
APH = 9;

% height
h                   = 10.^( -1:0.0001:(6-1e-6) )';

% obtain densities of atmospheric components
lat  = 0;       lon = 0;
currentYear = 2020;    dayOfYear = 1;     UTseconds = 0;
[~, rhos]           = atmosnrlmsise00(h,lat,lon,currentYear,dayOfYear,UTseconds,F107a,F107,APH);

% get density of air mixture up to 1e6 using nrlmsise m (then assume the
% density of space is 9e-27 kg/m^3
rho                 = [rhos(:,6); 9e-27; 9e-27];
h                   = [h; 1e6; 1e9];

% fit using a makima piecewise fit
system.atmosphericDensityMakimaPp  = makima( log10(h), log10(rho) );

if plotAirdensities
    
    hold on
    grid minor
    set(gca, 'YScale', 'log','XScale', 'log');
    plot( h, 10.^ppval(system.atmosphericDensityMakimaPp,log10(h)),'k-');
    plot( h([1:1000:end-3 end-2:end]), rho([1:1000:end-3 end-2:end]),'k*', ...
        MarkerSize=4 );
    
    xlabel("Height [$m$]",Interpreter='latex')
    ylabel("Air density [$kg/m^3$]",Interpreter='latex')
    set(gca, 'TickLabelInterpreter', 'latex');
    title("Earths atmospheric density as a function of height."+...
        newline+"Density of space is assumed to be $9e-27$ $kg/m^3$ at" + ...
        " height $>10^6$ $m$.",'Interpreter','latex')
    legend(["MAKIMA interpolation fit", "NRLMSISE-00 datapoints"], ...
        'Orientation','vertical','Location','best','Interpreter','latex')
    xlim(10.^[-2 inf])
    ylim([10^(-27), 1e1])
    f=gcf;
    set(f,'renderer','Painters')
    ax = gca;
    saveas(f,'NRLMSISE-00_Drag','epsc');
end
%%
% keep the workspace clean, delete what you don't need.
clear h rho rhos lat lon currentYear dayOfYear
clear UTseconds F107a F107 APH plotAirdensities;

% LOAD DATA ---------------------------------------------------------------
% selected systems
tblPropulsionSystems =  table2cell( readtable( "propulsionSystems.xlsx", "Sheet", "Data" ) );
% selected satellites
tblSatellites =         table2cell( readtable( "satellites.xlsx", "Sheet", "Data" ) );
% selected targets
tblTargets =            table2cell( readtable( "targets.xlsx", "Sheet", "Data" ) );

% Thrust = minimum( Tmax, Tmax * (Power available / Power required) )
%trueMaximumThrust = @( maxThrust, satellitePower, propulsionSystemPower )...
%    min( maxThrust , maxThrust * ( satellitePower / propulsionSystemPower ) );

% combine to scenarios
[nTargets,~] =      size( tblTargets );
[nSatellites,~] =   size( tblSatellites );
[nSystems,~] =      size( tblPropulsionSystems );

% SCENARIOS FORLOOPS TO WRITE TO SCENARIOS CELLTABLE
z = 1; tblScenarios = {};
for i = 1:nTargets; for j=1:nSatellites; for k=1:nSystems

    %SYSTEM DATA
    scenario.system                     = system;

    %TARGET DATA
    scenario.target.name                = tblTargets{i,1};
    scenario.target.state(1,1)          = conversion( tblTargets{i,2} ,"km","au");
    scenario.target.state(2,1)          = conversion( tblTargets{i,3} ,"km","au");
    scenario.target.state(3,1)          = conversion( tblTargets{i,4} ,"km","au");
    scenario.target.state(4,1)          = conversion( tblTargets{i,5} ,"km/s","au/y");
    scenario.target.state(5,1)          = conversion( tblTargets{i,6} ,"km/s","au/y");
    scenario.target.state(6,1)          = conversion( tblTargets{i,7} ,"km/s","au/y");
    scenario.target.approachDate        = datetime( string(tblTargets{i,8}) ,'Format', system.dateformat );
    
    %SATELLITE DATA
    scenario.satellite.mass             = tblSatellites{j,1};
    scenario.satellite.dryMass          = tblSatellites{j,7} + tblPropulsionSystems{k,5};
    scenario.satellite.P0               = tblSatellites{j,2};
    scenario.satellite.frontalArea      = conversion(tblSatellites{j,3} ,'m^2', 'au^2');
    scenario.satellite.dragCoefficient  = tblSatellites{j,4};
    scenario.satellite.apoapsisHeight   = conversion(tblSatellites{j,5},'km','au');
    scenario.satellite.apoapsis         = scenario.satellite.apoapsisHeight + system.bodies{1}.r;
    scenario.satellite.periapsisHeight  = conversion(tblSatellites{j,6},'km','au');
    scenario.satellite.periapsis        = scenario.satellite.periapsisHeight + system.bodies{1}.r;
    scenario.satellite.epsilon          = 1e-2;
    
    %PROPULSION SYSTEM DATA
    scenario.propulsionSystem.name      = tblPropulsionSystems{k,1};
    scenario.propulsionSystem.type      = tblPropulsionSystems{k,2};
    scenario.propulsionSystem.thrust    = conversion(tblPropulsionSystems{k,3}, 'm/s^2', 'au/y^2');
    scenario.propulsionSystem.Pprop     = tblPropulsionSystems{k,6};
    scenario.propulsionSystem.isp       = conversion(tblPropulsionSystems{k,4},'s','y');
    scenario.propulsionSystem.ve        = conversion(conversion(tblPropulsionSystems{k,4},'isp','effective exhaust velocity'),'m/s','au/y');
    
    %WRITE TO SCENARIO
    tblScenarios{z} = scenario;

    z = z + 1;
end;end;end
% keep the workspace clean, delete what you don't need.
clear i j k z scenario trueMaximumThrust;


[~,nScenarios] = size(tblScenarios);

% calculate Earth Escape for all scenarios following simple heuristic laws
numUniqueEarthEscapes = nScenarios/nTargets;

if ~isfile("results\results3D.mat")
    results = tblScenarios;
    save("results\results3D.mat","results");
end


for i = 1:numUniqueEarthEscapes
    load("results\results3D.mat","results");
    try
    if isfield(results{i},"earthEscape")
        if isfield(results{i}.earthEscape,"sol")
            continue;
        end
    end
    catch
    end
    % log message
    fprintf('Calculating Earth escape trajectory of the %s system on a %d kg satellite:\n',...
        tblScenarios{i}.propulsionSystem.name, ...
        tblScenarios{i}.satellite.mass );
    
    % create simulation model class object
    model  = EarthEscapeModel( tblScenarios{i} );

    % get initial state
    s0          = model.startingOrbit;
    
    % integrator ode settings
    odeSettings = odeset(  ...
            'RelTol', 5e-14, ...
            'AbsTol', 5e-14, ...
            'Events', @(t, s) model.odeEvents(t, s));
    
    % get solutions from ode89 integrator (RungeKutta 8 with global 
    % truncation error = O(h^9) due to variable step size)
    sol         = ode89(@(t, s) model.f(t, s), ...          % model
                        [ 0, max_mission_duration ], ...    %time interval
                        s0, ...                             % initial cond.
                        odeSettings);                       % ode settings
    

    % check solution found:
    fprintf('Result:\n')
    % mission fails if satellite does not escape Earth within max time
    if not(isfield(sol,'xe') && ~isempty(sol.xe))

        for j = 1:nTargets
            idsCases = i+(j-1)*nSystems*nSatellites;
            tblScenarios{idsCases}.earthEscape.sol = sol;
            tblScenarios{idsCases}.earthEscape.possible = false;
            results{idsCases} = tblScenarios{idsCases};
        end
        fprintf('Current scenario does not escape Earth within %d years!\n\n',max_mission_duration);

    % mission fails if satellite crashed into Earth, or runs out of
    % propellant:
    elseif ~any(sol.ie==1) % run out of fuel or crash into Earth 

        for j = 1:nTargets
            idsCases = i+(j-1)*nSystems*nSatellites;
            tblScenarios{idsCases}.earthEscape.sol = sol;
            tblScenarios{idsCases}.earthEscape.possible = false;
            results{idsCases} = tblScenarios{idsCases};
        end
        fprintf('Current scenario does not escape Earth!\n\n');

    % mission succes if Earth escape event occurs (same as else in this
    % case)
    else

        for j = 1:nTargets
            idsCases = i+(j-1)*nSystems*nSatellites;
            tblScenarios{idsCases}.earthEscape.sol = sol;
            tblScenarios{idsCases}.earthEscape.tf = sol.x(end);
            tblScenarios{idsCases}.earthEscape.sf = sol.y(:,end);
            tblScenarios{idsCases}.earthEscape.possible = true;
            results{idsCases} = tblScenarios{idsCases};
        end
        fprintf('Current scenario escapes Earth with %d kg mass left.\n\n',...
            sol.y(5,end)*tblScenarios{i}.satellite.mass)
    end
    save("results\results3D.mat","results");
end

% keep the workspace clean, delete what you don't need.
clear s0 sol simulation odeSettings;

%% Trajectory to near-Earth asteroid
for i=1:nScenarios
    % don't try to calculate the trajectory if Earth escape isn't possible.
    load("results\results3D.mat","results");
    try
    if isfield(results{i},"xopt")
        %continue;
    else
        tblScenarios{i} = results{i};
        if ~tblScenarios{i}.earthEscape.possible
            continue;
        end
    end
    catch
    end

    % load scenario to fly-by model
    model       = FlybyTrajectoryModel3D( tblScenarios{i} );

    % optimization settings
    nvars       = 7; % [trajectory duration and covars for: x, y, vx, vy, m]

    % optimization bounds for variables
    lb          = zeros( nvars, 1 ); ub = ones( nvars, 1 );
    
    % bounds for trajectory duration
    lb(1)       = .1;
    ub(1)       =  5;%max_mission_duration - tblScenarios{i}.earthEscape.tf;

    % Swarm size (higher sparm size results in better convergence to global
    % optim, but takes longer to calculate).
    swarmSize   = 2000 * nvars; % often enough, can be increased for better guess
    
    % particle swarm optimization settings
    psoOptions = optimoptions( ...
    'particleswarm',...
    'SwarmSize',                swarmSize,...
    'UseParallel',              true,... % Parallel proccessing TRUE
    'MaxIterations',            200,...  % A solution will likely not be found
    'Display',                  'iter',... % Show iterations
    'InertiaRange',             [.4 .9],... % Inertia settings from paper
    'SelfAdjustmentWeight',     1.49,... % Standard settings matlab
    'SocialAdjustmentWeight',   1.49,... % Standard settings matlab
    'MaxStallIterations',       20,...   % After 20 most likely cannot converge
    'FunctionTolerance',        conversion(1e5,'km','au')^2,... % (1e6 km in au)^2
    'ObjectiveLimit',           conversion(5e5,'km','au')^2,...
    'OutputFcn',@model.pswplotranges);   % (distance of 5e6 km in au)^2 is 
    % often sufficient for local gradient based optimization to converge.
    
    % local optimization using levenberg-marquandt method or damped least-squares
    fsolveopt = optimoptions("fsolve",...
            'Algorithm','levenberg-marquardt',...
            'UseParallel',true,... % parallel processing on
            'FiniteDifferenceType','forward',... % finite difference forward euler scheme sufficient
            'Display','iter',... % display iterations
            'StepTolerance',1e-20,... % step tolerance for variables 1e-20 can be necesary.
            'MaxFunctionEvaluations',1e5,... % converge not likely after 1e5 evaluations
            'MaxIterations',1e4,...  % converge not likely after 1e4 iterations
            'FunctionTolerance', conversion(1,'km','au'));%,... % 1 km tolerance
            %'ScaleProblem','Jacobian'); % use the pseudo-Jacobian to scale steps to take

    % calculate trajectory
    % global optimization using particle swarm: output --> guess
    [xopt,fval,exitflag,output]  = particleswarm( ...
      @(vars)model.OCP(vars, false, false), ...
      nvars, ...
      lb, ...
      ub, ...
      psoOptions);
    
    % use found solution from particle swarm as initial guess for local
    % optimization using levenberg-marquandt or damped least squarded
    % optimization:
    if fval<0.01 % if guess is good enough, else stop
        try
            [xopt2,fval,exitflag,output] = fsolve( ...
                @(vars)model.OCP(vars, true, false), ...
                xopt, ...
                fsolveopt);
        catch
            xopt2 = xopt;
        end
    else
        xopt2 = xopt;
    end
    % store found solution:
    tblScenarios{i}.xopt = xopt2;
    results{i} = tblScenarios{i};
    save("results\results3D.mat","results");
end
