% Use this file to run various scenarios and generate plots
%% set any common variables
concentration = 5; %grams/L influent
variables.Q = .47*1e-3;
variables.u =  25;
days =3;
variables.timepoints = 0:(days*60*24);%model x days of operation
season = "Uganda, Summer"; %This line labels the plots
csvfilename = "UGsumOutputsRevision.csv"; %File containing nitrogen data for each chamber
figfilename ="UGsumRevision"; %Figure name prefix
[tTemp,Tt]=TransferHeat([0,24*60]); %set the temperature function
variables.tTemp = tTemp;
variables.Tt = Tt;
%%
clc
close all

scenario = @(QoI,t) directgridScenario(QoI,t,variables.tTemp,variables.Tt)
scenName ="hybrid";

[t,err,removal,recovery]= getErrorBars(variables,scenario,chain,s2chain);
writematrix(err,scenName+csvfilename);
f=figure;
plotfixer
subplot(2,2,1)
hold on
trap = err(:,7:9);
plt= plotErrBars(t,err);
title("ECS operation")
ylabel("Relative N Concentration (% of A0)")
xlabel("Time (days)")
legend("A","C","T");


title("ECS operation")
ylabel("Relative N Concentration (% of A0)")
xlabel("time (days)")
legend("A","C","T");


%plot electricity 
subplot(2,2,2);
elements = length(t');
c_out = zeros(1,elements);
for i=1:elements
    c_out(i) = scenario("I",t(i));
end
plot(t/(60*24),c_out)
legend("ECS Input Current")
title("Electricity")
ylabel("mA/cm2")
xlabel("Time (days)")

%plot temperature profile
subplot(2,2,3);
plot(t/(60*24),scenario("T",t))
hold on
plot(t/(60*24),scenario("ST",t))
legend("Cathode Temp","Solar Panel Temp")
title("Temperatures")
ylabel("Kelvins")
xlabel("Time (days)")


subplot(2,2,4);
p_out = zeros(1,elements);
for i=1:elements
    p_out(i) = scenario("P",t(i));
end
gridcost = solarECSEconomics('W',p_out,t);

gramsNH3 = trap(:,2)*concentration; %Assume 5g/l influent
fertilizer = solarECSEconomics('NH3',gramsNH3,t);
fertilizerErr = solarECSEconomics('NH3',trap*concentration,t);
revenueErr = fertilizerErr-gridcost(2,:)';
plot(t/(60*24),gridcost(2,:),"DisplayName","Cumulative Electricty Cost")
hold on
fertRev = plotErrBars(t,fertilizerErr);
netRev = plotErrBars(t,revenueErr);
ylabel('$')
yyaxis right
plot(t/(60*24),gridcost(3,:),"DisplayName","Current Electricty cost", "Color",[0.4940 0.1840 0.5560])
ax=gca;
ax.YAxis(2).Color = [0.4940 0.1840 0.5560];
title("Economics")
legend("Cumulative Electricty Cost","Cumulative Fertlizer Revenue","Net Revenue","Current Electricty cost")
ylabel("$/kWh")
xlabel("Time (days)")

sgtitle('Direct PV + Grid Electricity - '+season) 
directgridNH3 = trap*concentration;
directgridrevenueErr= revenueErr;
directgridRemoval = removal;
directgridRecoveryDay =recovery*concentration/3;
figfile = figfilename+scenName;
saveas(f,figfile)

dollarGram = directgridrevenueErr(end,:)./(trap(end,:)*concentration)*1e3;
directgriddollarGram = [dollarGram(2), dollarGram(1)-dollarGram(3)];
out = {"Removal", removal;"RecoveryDay", recovery*concentration/3;"RevGram", [dollarGram(2), dollarGram(1)-dollarGram(3)] };


% mppt scenario
clc
%close all

scenario = @(QoI,t) mpptScenario(QoI,t,variables.tTemp,variables.Tt)
scenName = "mppt";
[t,err,removal,recovery]= getErrorBars(variables,scenario,chain,s2chain);
writematrix(err,scenName+csvfilename);
f=figure;
plotfixer
subplot(2,2,1)
hold on
trap = err(:,7:9);
plt= plotErrBars(t,err);
title("ECS operation")
ylabel("Relative N Concentration (% of A0)")
xlabel("Time (days)")
legend("A","C","T");


% close
% clc
%
subplot(2,2,2)
pow=[];
powcons =[];
for tp =t'
    P=scenario('P',tp);
    pow=[pow,P];
end
totalpow = cumtrapz(t,pow);
Vbat = 3.3471; %voltage for operation at 500 mA based on electrolysis curve (actually a bit lower than the data)
totalcons = cumtrapz(t,.5*Vbat*ones(size(t')));
plot(t/(60*24),totalpow,"DisplayName", "Energy Produced" )
hold on
plot(t/(60*24),totalcons, "DisplayName","Energy Used")
plot(t/(60*24),totalpow-totalcons,"DisplayName", "Energy Stored in Battery")
legend()
title("Electricity")
ylabel("Energy(J)")
xlabel("Time (days)")

%plot temperature profile

subplot(2,2,3);
plot(t/(60*24),scenario("T",t))
hold on
plot(t/(60*24),scenario("ST",t))
legend("Cathode Temp","Solar Panel Temp")
title("Temperatures")
ylabel("Kelvins")
xlabel("Time (days)")


subplot(2,2,4)
title("Economics")
ylabel("$")
xlabel("time (days)")
%assume electricity revenues are calculated as if we used and sold energy
%to the grid at the same time, revenue is the net gain

gridcost = solarECSEconomics('W',ones(size(t))*.5*Vbat,t);
gridgain = solarECSEconomics('sell',pow,t);

gramsNH3 = trap(:,2)*concentration; %Assume 5g/l influent
fertilizer = solarECSEconomics('NH3',gramsNH3,t);
fertilizerErr = solarECSEconomics('NH3',trap*concentration,t);
revenueErr = fertilizerErr-gridcost(2,:)'+gridgain(2,:)';
plot(t/(60*24),gridcost(2,:),"DisplayName","Cumulative Electricty Cost")
plot(t/(60*24),gridgain(2,:),"DisplayName","Cumulative Electricty Sale")
plot(t/(60*24),gridgain(2,:)-gridcost(2,:),"DisplayName","Net Electricty Sale")
hold on
fertRev = plotErrBars(t,fertilizerErr);
netRev = plotErrBars(t,revenueErr);
ylabel('$')
yyaxis right
plot(t/(60*24),gridcost(3,:),"DisplayName","Current Electricty cost", "Color",[0.4940 0.1840 0.5560])
plot(t/(60*24),gridgain(3,:),"DisplayName","Current Electricty price", "Color",[0.4940 0.1840 0.5560])
ax=gca;
ax.YAxis(2).Color = [0.4940 0.1840 0.5560];
title("Economics")
legend("Cumulative Electricty Cost","Cumulative Electricty Sale","Net Electricty Sale","Cumulative Fertlizer Revenue","Net Revenue","Current Electricty cost","Current Electricty price")
ylabel("$/kWh")
xlabel("Time (days)")

sgtitle('Controlled PV electricty - '+season) 
mpptNH3 = trap*concentration;
mpptrevenueErr= revenueErr;
mpptRemoval = removal;
mpptRecoveryDay =recovery*concentration/3;
dollarGram = mpptrevenueErr(end,:)./(trap(end,:)*concentration)*1e3;
mpptdollarGram = [dollarGram(2), dollarGram(1)-dollarGram(3)];
%out = {"Removal", removal;"RecoveryDay", recovery*concentration/3;"RevGram", [dollarGram(2), dollarGram(1)-dollarGram(3)] };
figfile = figfilename+scenName;
saveas(f,figfile)


% grid only scenario - 10 mA/cm^2
clc
%close all
idens = 10;
scenario = @(QoI,t) gridScenario(QoI,t,idens);
scenName="grid10";
f=figure;
plotfixer
subplot(2,2,1)
[t,err,removal,recovery]= getErrorBars(variables,scenario,chain,s2chain);
writematrix(err,scenName+csvfilename);
hold on
trap = err(:,7:9);
plt= plotErrBars(t,err);
title("ECS operation")
ylabel("Relative N Concentration (% of A0)")
xlabel("Time (days)")
legend("A","C","T");


%plot electricity 
subplot(2,2,2);
elements = length(t');
c_out = zeros(1,elements);
for i=1:elements
    c_out(i) = scenario("I",t(i));
end
plot(t/(60*24),c_out)
legend("ECS Input Current")
title("Electricity")
ylabel("mA/cm2")
xlabel("Time (days)")

%plot temperature profile
subplot(2,2,3);
plot(t/(60*24),scenario("T",t).*ones(size(t)))
legend("Cathode Temp")
title("Temperatures")
ylabel("Kelvins")
xlabel("Time (days)")


subplot(2,2,4);
p_out = zeros(1,elements);
for i=1:elements
    p_out(i) = scenario("P",t(i));
end
gridcost = solarECSEconomics('W',p_out,t);

gramsNH3 = trap(:,2)*concentration; %Assume 5g/l influent
fertilizer = solarECSEconomics('NH3',gramsNH3,t);
fertilizerErr = solarECSEconomics('NH3',trap*concentration,t);
revenueErr = fertilizerErr-gridcost(2,:)';
plot(t/(60*24),gridcost(2,:),"DisplayName","Cumulative Electricty Cost")
hold on
fertRev = plotErrBars(t,fertilizerErr);
netRev = plotErrBars(t,revenueErr);
ylabel('$')
yyaxis right
plot(t/(60*24),gridcost(3,:),"DisplayName","Current Electricty cost", "Color",[0.4940 0.1840 0.5560])
ax=gca;
ax.YAxis(2).Color = [0.4940 0.1840 0.5560];

title("Economics")
legend("Cumulative Electricty Cost","Cumulative Fertlizer Revenue","Net Revenue","Current Electricty cost")
ylabel("$/kWh")
xlabel("Time (days)")

sgtitle('Grid Electricity, i=10 mA/cm2 - '+season) 
grid10NH3 = trap*concentration;
grid10revenueErr= revenueErr;
grid10Removal = removal;
grid10RecoveryDay =recovery*concentration/3;
dollarGram = grid10revenueErr(end,:)./(trap(end,:)*concentration)*1e3;
grid10dollarGram = [dollarGram(2), dollarGram(1)-dollarGram(3)];
%out = {"Removal", removal;"RecoveryDay", recovery*concentration/3;"RevGram", [dollarGram(2), dollarGram(1)-dollarGram(3)] };
figfile = figfilename+scenName;
saveas(f,figfile);

%
idens = 15;
scenario = @(QoI,t) gridScenario(QoI,t,idens);
scenName="grid15";
f=figure;
plotfixer
subplot(2,2,1)
[t,err,removal,recovery]= getErrorBars(variables,scenario,chain,s2chain);
writematrix(err,scenName+csvfilename);
hold on
trap = err(:,7:9);
plt= plotErrBars(t,err);
title("ECS operation")
ylabel("Relative N Concentration (% of A0)")
xlabel("Time (days)")
legend("A","C","T");



%plot electricity 
subplot(2,2,2);
elements = length(t');
c_out = zeros(1,elements);
for i=1:elements
    c_out(i) = scenario("I",t(i));
end
plot(t/(60*24),c_out)
legend("ECS Input Current")
title("Electricity")
ylabel("mA/cm2")
xlabel("Time (days)")

%plot temperature profile
subplot(2,2,3);
plot(t/(60*24),scenario("T",t).*ones(size(t)))
legend("Cathode Temp")
title("Temperatures")
ylabel("Kelvins")
xlabel("Time (days)")


subplot(2,2,4);
p_out = zeros(1,elements);
for i=1:elements
    p_out(i) = scenario("P",t(i));
end
gridcost = solarECSEconomics('W',p_out,t);

gramsNH3 = trap(:,2)*concentration; %Assume 5g/l influent
fertilizer = solarECSEconomics('NH3',gramsNH3,t);
fertilizerErr = solarECSEconomics('NH3',trap*concentration,t);
revenueErr = fertilizerErr-gridcost(2,:)';
plot(t/(60*24),gridcost(2,:),"DisplayName","Cumulative Electricty Cost")
hold on
fertRev = plotErrBars(t,fertilizerErr);
netRev = plotErrBars(t,revenueErr);
ylabel('$')
yyaxis right
plot(t/(60*24),gridcost(3,:),"DisplayName","Current Electricty cost", "Color",[0.4940 0.1840 0.5560])
ax=gca;
ax.YAxis(2).Color = [0.4940 0.1840 0.5560];

title("Economics")
legend("Cumulative Electricty Cost","Cumulative Fertlizer Revenue","Net Revenue","Current Electricty cost")
ylabel("$/kWh")
xlabel("Time (days)")

sgtitle('Grid Electricity, i=15 mA/cm2 - '+season) 
grid15NH3 = trap*concentration;
grid15revenueErr= revenueErr;
grid15Removal = removal;
grid15RecoveryDay =recovery*concentration/3;
dollarGram = grid15revenueErr(end,:)./(trap(end,:)*concentration)*1e3;
grid15dollarGram = [dollarGram(2), dollarGram(1)-dollarGram(3)];
%out = {"Removal", removal;"RecoveryDay", recovery*concentration/3;"RevGram", [dollarGram(2), dollarGram(1)-dollarGram(3)] };
figfile = figfilename+scenName;
saveas(f,figfile);



%


% direct only scenario
clc
%close all


scenario = @(QoI,t) directScenario(QoI,t,variables.tTemp,variables.Tt)
scenName="direct";
f=figure;
plotfixer
subplot(2,2,1)
[t,err,removal,recovery]= getErrorBars(variables,scenario,chain,s2chain);
writematrix(err,scenName+csvfilename);
hold on
trap = err(:,7:9);
plt= plotErrBars(t,err);
title("ECS operation")
ylabel("Relative N Concentration (% of A0)")
xlabel("Time (days)")
legend("A","C","T");


%plot electricity 
subplot(2,2,2);
elements = length(t');
c_out = zeros(1,elements);
for i=1:elements
    c_out(i) = scenario("I",t(i));
end
plot(t/(60*24),c_out)
legend("ECS Input Current")
title("Electricity")
ylabel("mA/cm2")
xlabel("Time (days)")
@directScenario
%plot temperature profile
subplot(2,2,3);
plot(t/(60*24),scenario("T",t))
hold on
plot(t/(60*24),scenario("ST",t))
legend("Cathode Temp","Solar Panel Temp")
title("Temperatures")
ylabel("Kelvins")
xlabel("Time (days)")


subplot(2,2,4);
p_out = zeros(1,elements);
for i=1:elements
    p_out(i) = scenario("P",t(i));
end
gridcost = solarECSEconomics('W',p_out,t);

gramsNH3 = trap(:,2)*concentration; %Assume 5g/l influent
fertilizer = solarECSEconomics('NH3',gramsNH3,t);
fertilizerErr = solarECSEconomics('NH3',trap*concentration,t);
revenueErr = fertilizerErr-gridcost(2,:)';
plot(t/(60*24),gridcost(2,:),"DisplayName","Cumulative Electricty Cost")
hold on
fertRev = plotErrBars(t,fertilizerErr);
netRev = plotErrBars(t,revenueErr);
ylabel('$')
yyaxis right
plot(t/(60*24),gridcost(3,:),"DisplayName","Current Electricty cost", "Color",[0.4940 0.1840 0.5560])
ax=gca;
ax.YAxis(2).Color = [0.4940 0.1840 0.5560];

title("Economics")
legend("Cumulative Electricty Cost","Cumulative Fertlizer Revenue","Net Revenue","Current Electricty cost")
ylabel("$/kWh")
xlabel("Time (days)")

sgtitle('Direct PV Electricity - '+season) 
directNH3 = trap*concentration;
directrevenueErr= revenueErr;
directRemoval = removal;
directRecoveryDay =recovery*concentration/3;
dollarGram = directrevenueErr(end,:)./(trap(end,:)*concentration)*1e3;
directdollarGram = [dollarGram(2), dollarGram(1)-dollarGram(3)];
%out = {"Removal", removal;"RecoveryDay", recovery*concentration/3;"RevGram", [dollarGram(2), dollarGram(1)-dollarGram(3)] };
figfile = figfilename+scenName;
saveas(f,figfile)


aggregate = {"Metric", "Grid 10 mA", "Grid 15 mA", "Direct", "Hybrid", "MPPT";
    "Removal", grid10Removal,grid15Removal,directRemoval,directgridRemoval,mpptRemoval;
    "Recovery per Day", grid10RecoveryDay,grid15RecoveryDay,directRecoveryDay, directgridRecoveryDay,mpptRecoveryDay;
    "$/kg",grid10dollarGram,grid15dollarGram,directdollarGram, directgriddollarGram,mpptdollarGram}
writecell(aggregate,csvfilename)

%make economics plot
%close
%
f=figure;
scenName = "Allnoerror"

subplot(1,2,1)
sgtitle('Comparison of Operation Scenarios - '+season) 
plt = plot(t/(60*24),grid10NH3(:,2));
hold on
plt = plot(t/(60*24),grid15NH3(:,2));
plt = plot(t/(60*24),directNH3(:,2));
plt = plot(t/(60*24),directgridNH3(:,2));
plt = plot(t/(60*24),mpptNH3(:,2));

ylabel("Nitrogen Recovery (g N)")
xlabel("Time (days)")
title("Nitrogen Recovery")
legend("Grid electricty, i = 10 mA/cm2","Grid electricty, i = 15 mA/cm2", "Direct PV electricty","Combination PV and Grid electricty","Controlled PV electricty")
plotfixer


subplot(1,2,2)

g10 = plot(t/(60*24),grid10NH3(:,2));
hold on
g15 = plot(t/(60*24),grid15NH3(:,2));
direct = plot(t/(60*24),directNH3(:,2));
hybrid = plot(t/(60*24),directgridNH3(:,2));
mppt = plot(t/(60*24),mpptNH3(:,2));

plot(t/(60*24),zeros(size(t)),"k--",'HandleVisibility','off')
legend("Grid electricty, i = 10 mA/cm2","Grid electricty, i = 15 mA/cm2", "Direct PV electricty","Combination PV and Grid electricty","Controlled PV electricty")
plotfixer
title("Revenues")
ylabel("Net Revenue ($)")
xlabel("Time (days)")
figfile = figfilename+scenName;
saveas(f,figfile);

%% cleaner figure?

f=figure;
scenName = "All"

subplot(1,2,1)

sgtitle('Comparison of Operation Scenarios - '+season) 

g10 = plot(t,grid10NH3(:,2));
hold on
g15 = plot(t,grid15NH3(:,2));
direct = plot(t,directNH3(:,2));
hybrid = plot(t,directgridNH3(:,2));
mppt = plot(t,mpptNH3(:,2));

g10up = plot(t,grid10NH3(:,1),"-.");
g15up = plot(t,grid15NH3(:,1),"-.");
directup = plot(t,directNH3(:,1),"-.");
hybridup = plot(t,directgridNH3(:,1),"-.");
mpptup = plot(t,mpptNH3(:,1),"-.");

g10down = plot(t,grid10NH3(:,3),"-.");
g15down = plot(t,grid15NH3(:,3),"-.");
directdown = plot(t,directNH3(:,3),"-.");
hybriddown = plot(t,directgridNH3(:,3),"-.");
mpptdown = plot(t,mpptNH3(:,3),"-.");

g10up.SeriesIndex = g10.SeriesIndex;
g10down.SeriesIndex = g10.SeriesIndex;
g15up.SeriesIndex = g15.SeriesIndex;
g15down.SeriesIndex = g15.SeriesIndex;
directup.SeriesIndex = direct.SeriesIndex;
directdown.SeriesIndex = direct.SeriesIndex;
hybridup.SeriesIndex = hybrid.SeriesIndex;
hybriddown.SeriesIndex = hybrid.SeriesIndex;
mpptup.SeriesIndex = mppt.SeriesIndex;
mpptdown.SeriesIndex = mppt.SeriesIndex;


ylabel("grams N recovered")
xlabel("time (days)")
title("Nitrogen Recovery")
legend()
plotfixer


subplot(1,2,2)

title("Revenues")
g10 = plot(t,grid10revenueErr(:,2),'k');
hold on
g15 = plot(t,grid15revenueErr(:,2),'k');
direct = plot(t,directrevenueErr(:,2),'k');
hybrid = plot(t,directgridrevenueErr(:,2),'k');
mppt = plot(t,mpptrevenueErr(:,2),'k');

g10up = plot(t(1:100:end),grid10revenueErr(1:100:end,1),"--");
g15up = plot(t(1:100:end),grid15revenueErr(1:100:end,1),"--");
directup = plot(t(1:100:end),directrevenueErr(1:100:end,1),"--");
hybridup = plot(t(1:100:end),directgridrevenueErr(1:100:end,1),"--");
mpptup = plot(t(1:100:end),mpptrevenueErr(1:100:end,1),"--");

g10down = plot(t(1:100:end),grid10revenueErr(1:100:end,3),"--");
g15down = plot(t(1:100:end),grid15revenueErr(1:100:end,3),"--");
directdown = plot(t(1:100:end),directrevenueErr(1:100:end,3),"--");
hybriddown = plot(t(1:100:end),directgridrevenueErr(1:100:end,3),"--");
mpptdown = plot(t(1:100:end),mpptrevenueErr(1:100:end,3),"--");

g10up.SeriesIndex = g10.SeriesIndex;
g10down.SeriesIndex = g10.SeriesIndex;
g15up.SeriesIndex = g15.SeriesIndex;
g15down.SeriesIndex = g15.SeriesIndex;
directup.SeriesIndex = direct.SeriesIndex;
directdown.SeriesIndex = direct.SeriesIndex;
hybridup.SeriesIndex = hybrid.SeriesIndex;
hybriddown.SeriesIndex = hybrid.SeriesIndex;
mpptup.SeriesIndex = mppt.SeriesIndex;
mpptdown.SeriesIndex = mppt.SeriesIndex;

yyaxis right
ax=gca;
ax.YAxis(2).Color = [0.4940 0.1840 0.5560];



plot(t/(60*24),zeros(size(t)),"k--",'HandleVisibility','off')
legend("Grid electricty, i = 10 mA/cm2","Grid electricty, i = 15 mA/cm2", "Direct PV electricty","Combination PV and Grid electricty","Controlled PV electricty")
plotfixer
ylabel("Net Revenue ($)")
xlabel("time (days)")
figfile = figfilename+scenName;

%%
function [t,err,removal,recovery]= getErrorBars(variables,scenario,chain,s2chain)
%currently only returns error on last time point
    progress = waitbar(1, 'hello')
chain = chain(1:100:end,:);
s2chain = s2chain(1:100:end); %fewer realizations
nt = length(variables.timepoints);

realizationsA=zeros(nt,size(chain,1));
    realizationsC=realizationsA;
    realizationsT=realizationsA;
    realizationsTAN=realizationsA;

    for j =1:size(chain,1)
        msg = "Percent complete" + string(j/size(chain,1)*100)+"%";
        waitbar(j/size(chain,1),progress,msg)
        c= chain(j,:);
        paramIN = c;
        [t,A_t] = extrapolateNumericalInt(paramIN,variables,scenario);

        s2 =sqrt(s2chain(j));
        realizationsA(:,j) = A_t(:,1)+randn(nt,1)*s2;
        realizationsC(:,j)= A_t(:,2)+randn(nt,1)*s2;
        realizationsT(:,j)= A_t(:,3)+randn(nt,1)*s2;
        realizationsTAN(:,j)=realizationsA(j)+realizationsC(j)+realizationsT(j);

        % realizationsA(j) = A_t(end,1)+randn(1,1)*s2;
        % realizationsC(j)= A_t(end,2)+randn(1,1)*s2;
        % realizationsT(j)= A_t(end,3)+randn(1,1)*s2;
        % realizationsTAN(j)=realizationsA(j)+realizationsC(j)+realizationsT(j);
    end
    intervalA = prctile(realizationsA,[95,50,5],2);
    intervalC = prctile(realizationsC,[95,50,5],2);
    intervalT = prctile(realizationsT,[95,50,5],2);
    intervalTAN= prctile(realizationsTAN,[95,50,5],2);
%err= [intervalA(1),intervalA(3);intervalC(1),intervalC(3);intervalT(1),intervalT(3);intervalTAN(1),intervalTAN(3)];
err= [intervalA, intervalC,intervalT];
removal = [(1-intervalA(end,2)),(intervalA(end,1)-intervalA(end,3))]*100;
recovery = [intervalT(end,2),(intervalT(end,1)-intervalT(end,3))];
close(progress)
end

function plot = plotErrBars(t,err)
thick = 1.5;
hold on
colorlist = {"b", [0.6350 0.0780 0.1840],[0.4660 0.6740 0.1880],"k"};
nLines = floor(size(err,2)/3);
for i=1:nLines
    if nLines > 1
        avg =(i*3-1);
        up = (i*3-2);
        down = i*3;
    plot = shadedErrorBar(t/(60*24),err(:,avg)',[abs(err(:,up)-err(:,avg))';abs(err(:,avg)-err(:,down))'],'lineProps',{"Color",colorlist{i},'LineWidth',thick});
    else
        plot = shadedErrorBar(t/(60*24),err(:,i+1)',[abs(err(:,i)-err(:,i+1))';abs(err(:,i+1)-err(:,i+2))'],'lineProps',{"LineWidth",.5});
    end
end
end
