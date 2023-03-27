%{
This script creates the figures for a low thrust trajectory using 
Global Optimization Toolbox's Partical Swarm Optimization
while using an optimal control strategy using the calculus of variations
indirect single shooting method

dependencies:
    Aerospace Toolbox
    Global Optimization Toolbox
    Parallel Processing Toolbox
    MICE
    JPL Horizons
    SPICE KERNELS

Parts:
    Cases: defines the different scenarios and combinations that should be
    simulated.

    Constants: defines the constants of the trajectory such as bodies,
    target, spacecraft and the trajectories of bodies and target.
    
    Then the initial orbit is defined as 90000 - 295 km altitude SSGTO
    around earth in the EarthEscapeModel.m class.

    Finally, the states of the satellite are defined as pos/vel/mass and
    their corresponding lagrange multipliers as used in the Hamiltonian.

    FlybyTrajectoryModel.m class handles the definition of all the 
    dynamics, control and optimization cost function.
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

% get number of unique targets
load("results\results3D.mat","results");
nScenarios = length(results);
targets={};for i = 1:nScenarios; targets{i} = results{i}.target.name;end
targets=unique(targets,'stable');
nTargets=length(targets);
systems={};for i = 1:nScenarios/nTargets; systems{i} = results{i}.propulsionSystem.name;end; systems=unique(systems,'stable');
nSystems=length(systems);

for i=1:nSystems
    for j=(1:nTargets)-1
        resultsPerSystem{i,j+1} = results{i+j*9};
    end
end

%% Calculate trajectories
for iSys = 1:nSystems
    for j = 1:nTargets
        % create simulation model
        simulation = FlybyTrajectoryModel3D( resultsPerSystem{ iSys, j } );
        
        % calculate interplanetary trajectory
        [fval, t, y, alpha, u] = simulation.OCP( resultsPerSystem{ iSys, j }.xopt, true, true );
        
        % write results in km and km/s
        yb =            zeros(size(y));
        yb(1:3,:) =     conversion(y(1:3,:), 'au', 'km');
        yb(4:6,:) =     conversion(y(4:6,:), 'au/y', 'km/s');
        yb(7:14,:) =    y(7:14,:);
        
        % write to resultsPerSystem
        resultsPerSystem{ iSys, j }.trajectoryTStart =  simulation.tStart; 
        resultsPerSystem{ iSys, j }.trajectoryTime =    t;
        resultsPerSystem{ iSys, j }.trajectory =        y;
        resultsPerSystem{ iSys, j }.trajectorykm =      yb;
        resultsPerSystem{ iSys, j }.trajectoryControl = u.*alpha;

        cspice_furnsh( resultsPerSystem{ iSys, j }.system.kernels )
        resultsPerSystem{ iSys, j }.system.sol_ets        = zeros(1,length(t));
        for n=1:length(t)
            resultsPerSystem{ iSys, j }.system.sol_ets(n) = cspice_str2et(char(resultsPerSystem{ iSys, j }.trajectoryTStart + years(t(n))));
        end
        [resultsPerSystem{ iSys, j }.earth.sol_states, ~] = cspice_spkezr(...
        resultsPerSystem{ iSys, j }.system.bodies{1}.name,...
        resultsPerSystem{ iSys, j }.system.sol_ets,...
        resultsPerSystem{ iSys, j }.system.frame,...
        resultsPerSystem{ iSys, j }.system.abcorr,...
        resultsPerSystem{ iSys, j }.system.observer);
        cspice_kclear;
    end
end
%% create figures x,y (z?)
% this really depends on the targets change for other targets!!
symbolList = ["v","square","diamond","o"];
colorList = ["black","black","black","black"];
names = ["2020QN1", "163693", "2017WV13", "2012BX34"];


for range = 1:nTargets
    j = range;
for iSys = 1:nSystems
figure;
t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile;
hold on;

% Sun
sunPlot = plot3(0,0,0,'o','Color',[0.9290 0.6940 0.1250],'MarkerFaceColor',[0.9290 0.6940 0.1250],'MarkerSize',10,'DisplayName','Sun');

% trajectory
y = resultsPerSystem{ iSys, j }.trajectory;

% thrusting action
u = resultsPerSystem{ iSys, j }.trajectoryControl;

% Earth
es =            resultsPerSystem{ iSys, j }.earth.sol_states;

% Earth in AU
esAU =          zeros( size( es ) );
esAU(1:3,:) =   conversion(es(1:3,:),'km', 'au');
esAU(4:6,:) =   conversion(es(4:6,:),'km/s', 'au/year');

% plot Earth
plot3(esAU(1,1), esAU(2,1),esAU(3,1),'Marker','x', 'LineStyle','none','Color',[0 0.4470 0.7410],'DisplayName', 'Earth');
plot3(   esAU(1,:), esAU(2,:),esAU(3,:),'--','Color',[0 0.4470 0.7410],'LineWidth',.5,'DisplayName', 'Earth');
earthPlot = plot3(esAU(1,end), esAU(2,end),esAU(3,end),'Marker','o', 'LineStyle','none','Color',[0 0.4470 0.7410],'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName', 'Earth');

% plot trajectory spacecraft
plot3(resultsPerSystem{ iSys, j }.target.state(1,:),...
     resultsPerSystem{ iSys, j }.target.state(2,:),...
     resultsPerSystem{ iSys, j }.target.state(3,:),symbolList(j),'Color','black','LineWidth',1);

% plot thicker lines when propulsion system is turned on
for b = 1:(length(y(1,:))-1)
    if norm(u(1:3,b))>0.5
         qw{j} = plot3(y(1,b:b+1),y(2,b:b+1),y(3,b:b+1),'-','Color',colorList(j),'LineWidth',2,'DisplayName', names(j));
    else
         plot3(y(1,b:b+1),y(2,b:b+1),y(3,b:b+1),'-','Color',colorList(j),'LineWidth',.5,'DisplayName', names(j));
    end
end

% legend
legend([sunPlot,earthPlot,qw{range}],'location','southoutside','Orientation','horizontal');
axh = gca(); % Get handle to axis
daspect([1 1 0.05])
pbaspect([1 1 1])
grid(axh,'minor')
xlim(axh,xlim*1.05)    
ylim(axh,ylim*1.05)
xlabel('$X [AU]$','Interpreter','latex')
ylabel('$Y [AU]$','Interpreter','latex')
zlabel('$Z [AU]$','Interpreter','latex')
title({"Cruise from Earth to asteroid using ", resultsPerSystem{ iSys, j }.propulsionSystem.name},'Interpreter','latex')
view(-45,35.264);

f=gcf;
set(f,'renderer','Painters')
ax = gca;
set(gca,'TickLabelInterpreter','latex');
rangeString = strjoin(string(range),'');
saveas(f,strjoin([resultsPerSystem{ iSys, j }.propulsionSystem.name, rangeString],''),'epsc');
end
end
%%
for iSys = 1:nSystems
figure;
t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
nexttile;
hold on;
earthPlot = plot3(0, 0, 0,'Marker','o', 'LineStyle','none','Color',[0 0.4470 0.7410],'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName', 'Earth');
y = resultsPerSystem{ iSys, j }.earthEscape.sol.y;
yS = [];
yS(1:3,:) = conversion([y(1,:) .* cos(y(2,:)); y(1,:) .* sin(y(2,:)); y(1,:) .* 0],'au','km');

qw = plot3(yS(1,:),yS(2,:),yS(3,:),'-','Color',"black",'LineWidth',.5,'DisplayName','Earth escape trajectory');

legend([earthPlot,qw],'location','southoutside','Orientation','horizontal');
axh = gca(); % Get handle to axis
%pbaspect([1 1 .2])
axis equal;
grid(axh,'minor')
%xlim(axh,xlim*1.05)    
%ylim(axh,ylim*1.05)
daspect([1 1 1])
pbaspect([1 1 .2])
xlabel('$X$ (km)','Interpreter','latex')
ylabel('$Y$ (km)','Interpreter','latex')
zlabel('$Z$ (km)','Interpreter','latex')
title({"Earth escape trajectory using ", resultsPerSystem{ iSys, j }.propulsionSystem.name},'Interpreter','latex')
view(2);
nexttile;
hold on;
earthPlot = plot3(0, 0, 0,'Marker','o', 'LineStyle','none','Color',[0 0.4470 0.7410],'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName', 'Earth');
qw = plot3(yS(1,:),yS(2,:),yS(3,:),'-','Color',"black",'LineWidth',.1,'DisplayName','Earth escape trajectory');
axh = gca(); % Get handle to axis
axis equal
grid(axh,'minor')
ax = gca;
ax.XLim = 1e5*[-1, 1] - 45000;
ax.YLim = 1e5*[-1, 1];
daspect([1 1 1])
pbaspect([1 1 .2])
xlabel('$X$ (km)','Interpreter','latex')
ylabel('$Y$ (km)','Interpreter','latex')
zlabel('$Z$ (km)','Interpreter','latex')
view(2);
f=gcf;
set(f,'renderer','Painters')
ax = gca;
set(gca,'TickLabelInterpreter','latex');
rangeString ="EarthEscape";
saveas(f,strjoin([resultsPerSystem{ iSys, j }.propulsionSystem.name, rangeString],''),'epsc');
end

%% create table data
targ = 1;

columnNames = ["System", "startDate", "ToF (days)", "Mp (kg)", "Mdry (kg)", "Distance (km)"];
clear colSystems colStartDate colToF colMp colMdry colD
for i = 1:nSystems
colSystems{i,1} =     resultsPerSystem{i,targ}.propulsionSystem.name;
colMdry{i,1} =        resultsPerSystem{i,targ}.satellite.dryMass - 3.4390;

targ=1;
colStartDate1{i,1} =   resultsPerSystem{i,targ}.trajectoryTStart - years(resultsPerSystem{i,targ}.earthEscape.tf);
colToF1{i,1} =         days(years(resultsPerSystem{i,targ}.xopt(1) + resultsPerSystem{i,targ}.earthEscape.tf));
colMp1{i,1} =          resultsPerSystem{i,targ}.satellite.mass * (1 - resultsPerSystem{i,targ}.trajectory(7,end));
colD1{i,1} =           conversion(norm(resultsPerSystem{i,targ}.trajectory(1:3,end) - resultsPerSystem{i,targ}.target.state(1:3,end)),'au','km');
colM1{i,1} =          resultsPerSystem{i,targ}.satellite.dryMass - 3.4390 + resultsPerSystem{i,targ}.satellite.mass * (1 - resultsPerSystem{i,targ}.trajectory(7,end));

targ=2;
colStartDate2{i,1} =   resultsPerSystem{i,targ}.trajectoryTStart - years(resultsPerSystem{i,targ}.earthEscape.tf);
colToF2{i,1} =         days(years(resultsPerSystem{i,targ}.xopt(1) + resultsPerSystem{i,targ}.earthEscape.tf));
colMp2{i,1} =          resultsPerSystem{i,targ}.satellite.mass * (1 - resultsPerSystem{i,targ}.trajectory(7,end));
colD2{i,1} =           conversion(norm(resultsPerSystem{i,targ}.trajectory(1:3,end) - resultsPerSystem{i,targ}.target.state(1:3,end)),'au','km');
colM2{i,1} =          resultsPerSystem{i,targ}.satellite.dryMass - 3.4390 + resultsPerSystem{i,targ}.satellite.mass * (1 - resultsPerSystem{i,targ}.trajectory(7,end));

targ=3;
colStartDate3{i,1} =   resultsPerSystem{i,targ}.trajectoryTStart - years(resultsPerSystem{i,targ}.earthEscape.tf);
colToF3{i,1} =         days(years(resultsPerSystem{i,targ}.xopt(1) + resultsPerSystem{i,targ}.earthEscape.tf));
colMp3{i,1} =          resultsPerSystem{i,targ}.satellite.mass * (1 - resultsPerSystem{i,targ}.trajectory(7,end));
colD3{i,1} =           conversion(norm(resultsPerSystem{i,targ}.trajectory(1:3,end) - resultsPerSystem{i,targ}.target.state(1:3,end)),'au','km');
colM3{i,1} =          resultsPerSystem{i,targ}.satellite.dryMass - 3.4390 + resultsPerSystem{i,targ}.satellite.mass * (1 - resultsPerSystem{i,targ}.trajectory(7,end));

targ=4;
colStartDate4{i,1} =   resultsPerSystem{i,targ}.trajectoryTStart - years(resultsPerSystem{i,targ}.earthEscape.tf);
colToF4{i,1} =         days(years(resultsPerSystem{i,targ}.xopt(1) + resultsPerSystem{i,targ}.earthEscape.tf));
colMp4{i,1} =          resultsPerSystem{i,targ}.satellite.mass * (1 - resultsPerSystem{i,targ}.trajectory(7,end));
colD4{i,1} =           conversion(norm(resultsPerSystem{i,targ}.trajectory(1:3,end) - resultsPerSystem{i,targ}.target.state(1:3,end)),'au','km');
colM4{i,1} =          resultsPerSystem{i,targ}.satellite.dryMass - 3.4390 + resultsPerSystem{i,targ}.satellite.mass * (1 - resultsPerSystem{i,targ}.trajectory(7,end));

end
tblData1 = table(colSystems,colStartDate1,colToF1,colMp1,colD1,colM1,...
                            colStartDate2,colToF2,colMp2,colD2,colM2,...
                            colStartDate3,colToF3,colMp3,colD3,colM3,...
                            colStartDate4,colToF4,colMp4,colD4,colM4,...
                            colMdry);
%%
f = figure('Renderer', 'painters', 'Position', [10 10 700 500]);
t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
nexttile;
title('Time of flight (ToF) and Total mass (kg) for each propulsion system to fly by selected asteroid.',Interpreter='latex');
hold on;
bar(reordercats(categorical(tblData1.colSystems),tblData1.colSystems),years(days([tblData1.colToF1{:};tblData1.colToF2{:};tblData1.colToF3{:};tblData1.colToF4{:}])))
yline(5,'--');
%leg = legend(["2020QN1", "163693", "2017WV13", "2012BX34",""],"Location","bestoutside",Orientation="horizontal",Interpreter="latex");
%title(leg,"Asteroid fly by");
grid minor;
ax = gca;
set(gca,'XTick',[]);
set(gca,'YTick',0:1:10);
yAx = get(gca,'YAxis');
yl = ylim;
yAx.MinorTickValues=yl(1):0.5:yl(2);
ylabel('ToF (years)','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');

%t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile;
%title('Total mass for each propulsion system per asteroid.',Interpreter='latex');
hold on;
bar(reordercats(categorical(tblData1.colSystems), tblData1.colSystems),[tblData1.colM1{:};tblData1.colM2{:};tblData1.colM3{:};tblData1.colM4{:}],0.8,"grouped")
bar(reordercats(categorical(tblData1.colSystems), tblData1.colSystems),[tblData1.colMdry{:};tblData1.colMdry{:};tblData1.colMdry{:};tblData1.colMdry{:}],0.8,"grouped","black","FaceAlpha",0.4,"EdgeColor","none")

leg = legend(["Asteroid 2020QN1", "Asteroid 163693", "Asteroid 2017WV13", "Asteroid 2012BX34","$M_{dry}$","","",""],"Location","southoutside",Orientation="horizontal",Interpreter="latex");
grid minor;
ax = gca;

ylabel('$M_{dry}+M_p$ (kg)','Interpreter','latex');
set(gca,'YTick',0:2:14);
yAx = get(gca,'YAxis');
yl = ylim;
yAx.MinorTickValues=yl(1):0.5:yl(2);
set(gca,'TickLabelInterpreter','latex');
exportgraphics(f, 'LowThrustSystemsSummary.pdf');