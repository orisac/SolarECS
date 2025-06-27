% Required package: The DRAM package is available at the websites
%           https://wiki.helsinki.fi/display/inverse/Adaptive+MCMC
%           http://helios.fmi.fi/~lainema/mcmc/

close all; clc

% Load file containing experimental data for each experimental conditions,
% including mean values and error bars
load("Data/ECSData.mat") 

% Normalize all concentrations with respect to the initial anode
% concentration
errmA5=mA5/mA5(1,1).*sqrt((errmA5(1,1)/mA5(1,1))^2+(errmA5./mA5).^2);
mA5=mA5/mA5(1,1);
errmA10=mA10/mA10(1,1).*sqrt((errmA10(1,1)/mA10(1,1))^2+(errmA10./mA10).^2);
mA10=mA10/mA10(1,1);
errmA15=mA15/mA15(1,1).*sqrt((errmA15(1,1)/mA15(1,1))^2+(errmA15./mA15).^2);
mA15=mA15/mA15(1,1);
errmA20=mA20/mA20(1,1).*sqrt((errmA20(1,1)/mA20(1,1))^2+(errmA20./mA20).^2);
mA20=mA20/mA20(1,1);
errT30=T30/T30(1,1).*sqrt((errT30(1,1)/T30(1,1))^2+(errT30./T30).^2);
T30=T30/T30(1,1);
errT40=T40/T40(1,1).*sqrt((errT40(1,1)/T40(1,1))^2+(errT40./T40).^2);
T40=T40/T40(1,1);
errT50=T50/T50(1,1).*sqrt((errT50(1,1)/T50(1,1))^2+(errT50./T50).^2);
T50=T50/T50(1,1);
errT60=T60/T60(1,1).*sqrt((errT60(1,1)/T60(1,1))^2+(errT60./T60).^2);
T60=T60/T60(1,1);



%Package data in accessible formats
expData =[mA5;mA10;mA15;mA20;T30;T40;T50;T60];
TAN = sum(expData,2); %Mass Balance
TANData = {sum(mA5,2),sum(mA10,2),sum(mA15,2),sum(mA20,2),sum(T30,2),sum(T40,2),sum(T50,2),sum(T60,2)};
errorData = [errmA5;errmA10;errmA15;errmA20;errT30;errT40;errT50;errT60];
errTAN = sqrt(sum(errorData.^2,2));
TANerr = {sqrt(sum(errmA5.^2,2)),sqrt(sum(errmA10.^2,2)),sqrt(sum(errmA15.^2,2)),sqrt(sum(errmA20.^2,2)),sqrt(sum(errT30.^2,2)),sqrt(sum(errT40.^2,2)),sqrt(sum(errT50.^2,2)),sqrt(sum(errT60.^2,2))};
dataC = {mA5,mA10,mA15,mA20,T30,T40,T50,T60};
dataErr = {errmA5,errmA10,errmA15,errmA20,errT30,errT40,errT50,errT60};
pH =[mA5ph,mA10ph,mA15ph,mA20ph,T30ph,T40ph,T50ph,T60ph];

nC = size(conditions,1);

data = {};
data.xdata = timepoints;
data.udata = expData;
data.error = errorData;
n=numel(expData);
data.conditions = conditions;
data.ph = pH;
%% Construct guesses and sensitivites
p=5; %number of parameters


param0=[285.64, -4005.6, 4842.86, -4100, 13.47]; %Educated initial guess derived from 1st principles


data.fix =0;
data.errScale =0;
obj = @(params) sumSquaresNumericalInt(params, data);
options = optimset('MaxIter',9000*p,'MaxFunEvals', 9000*p,'TolX',1e-8,'Display', 'notify');%5x the default max iterations

%Improve initial guess for DRAM search using computationally inexpensive
%fminsearch
[paramsI,fval,exitflag] = fminsearch(obj, param0)%,options);
prun = paramsI;
for i = 1:15
    param0 = paramsI;
    flast = fval;
    [paramsI,fval,exitflag] = fminsearch(obj, param0,options)
    if abs(fval-flast)<1e-4 %run until convergence 
        break
    end
end


% Construct sensitivities to create covariance matrix
close all
sensfig =figure;
c0=1;
modelA=zeros(8*nC,p);
modelC=zeros(8*nC,p);
modelT=zeros(8*nC,p);
nC = size(conditions,1);
variables.timepoints = timepoints;
senstime = timepoints(2:end); %dont include initial val since there is no sensitivity
ss_sensitivity = zeros(1,8);
variables.fix = 0;
data.errScale =0;
for j = 1:p
    sens_param = paramsI;
    sens_param(j) = 1.1*paramsI(j);
    model_plus = [];
    model_plusSum=0;
    for i = 1:nC
        variables.I = conditions(i,1);
        variables.T = conditions(i,2);
        variables.u = conditions(i,3);
        A_t = plotNumericalInt(sens_param,variables); %use secondary function for more accessible outputs
        model_plusSum = model_plusSum+A_t(2:9,:);
        model_plus = [model_plus;A_t(2:9,:)];
        ss_plus = sumSquaresNumericalInt(sens_param, data);
    end
    model_minus = [];
    model_minusSum =0;
    sens_param = paramsI;
    sens_param(j) = .9*paramsI(j);
    for i = 1:nC
        variables.I = conditions(i,1);
        variables.T = conditions(i,2);
        variables.u = conditions(i,3);
        A_t = plotNumericalInt(sens_param,variables);
        model_minusSum = model_minusSum+A_t(2:9,:);
        model_minus = [model_minus;A_t(2:9,:)];
        ss_minus = sumSquaresNumericalInt(sens_param, data);
    end
    dMda = (model_plus-model_minus)/ .2*paramsI(j);
    dMdaPlot = (model_plusSum-model_minusSum)/ .2*paramsI(j);
    ss_sensitivity(j) = (ss_plus-ss_minus)/ .2*paramsI(j);
    modelA(:,j)= dMda(:,1);
    modelC(:,j)= dMda(:,2);
    modelT(:,j)= dMda(:,3);
    disp = 'a'+string(j);
    %plot local sensitivity analysis
    for chamber =1:3
    sensfig(chamber)= subplot(1,3,chamber);
    set(sensfig(chamber),'colororder',turbo(p))
    hold on
    plot(timepoints(2:9),dMdaPlot(:,chamber),"DisplayName",disp)
    end
    
end
plotfixer
sens_matA = modelA'*modelA;
sens_matC = modelC'*modelC;
sens_matT = modelT'*modelT;
ss_sensitivity; %Remove colon to see sensitivity summary
% residuals
modelR =[];
variables.timepoints = timepoints;
for i = 1:size(conditions,1)
    variables.I = conditions(i,1);
    variables.T = conditions(i,2);
    variables.u = conditions(i,3);
modelR =[modelR; plotNumericalInt(paramsI,variables)];
end
 res = modelR - expData;
 resT = res(:,3);%sum(res,1)
% Construct the covariance matrix

sigma2 = (1/(n-p))*(resT'*resT);
V = sigma2*pinv(sens_matA*sens_matA');
D= diag(abs(paramsI));

%% DRAM

% Set the options employed in DRAM.
clear data model options
covMat =sigma2*sens_matA;
data.xdata =timepoints ;
data.udata = expData;
data.conditions = conditions;
data.error = errorData;
tmin = paramsI;
fixed = 0;
data.fix = fixed;
data.errScale =0;

tcov = covMat;

paramsN = { % parameter name, initial value and minimum
           {'a1', tmin(1), 0, inf}
           {'a2', tmin(2), -inf, 0}
           {'a3', tmin(3), 0, inf}
           {'a4', tmin(4), -inf, 0}
           {'a5', tmin(5), 0, inf}};


model.ssfun   = @sumSquaresNumericalInt;    % function that computes sum of squares error
model.sigma2  = sigma2;    % sigma2 of jumping distribution
options.nsimu = 1e5;    % number of steps
options.qcov  = covMat;%results.qcov;%D;    % covariance matrix of jumping distribution
options.updatesigma = 1;
model.N       = n;    % number of data points

% Run DRAM to construct the chains (stored in chain) and measurement
% variance (stored in s2chain).
[results, chain, s2chain] = mcmcrun(model, data, paramsN, options);


% Plot
close all
figure(1); clf
mcmcplot(chain, [],results,'chainpanel');

figure(2); clf
mcmcplot(chain, [],results,'pairs');

figure(3); clf
mcmcplot(chain, [],results,'denspanel');

cov(chain)
chainstats(chain,results)


%% plot predicition intervals


%close all
degC = [char(176), 'C'];
variables.pfixed = paramsI;%(pnon);
variables.timepoints = timepoints;
timepoint = timepoints/60; %convert minutes to hours for plotting
variables.fix = fixed;
hf=figure;
fontsize(scale=1.2)
meanFit = [];
thick = 1.5;
for i = 1:size(conditions,1)

    if i == 5 %plot current and temperature experiments in seperate subplots
        legend("Anode data", "Cathode data", "Trap Data", "Total Data", "Anode Model", "Cathode Model", "Trap Model", "Total Model")
        sgtitle("Fitting Set")
        hf=figure;
        hold on
        fontsize(scale=1.2)
    end
    if i >4
        hf(i-4)= subplot(2,2,i-4);
    else
        hf(i)= subplot(2,2,i);
    end
    %hf(i)= subplot(2,2,i);
   
    g=errorbar(timepoint,dataC{i},dataErr{i},"o",'LineWidth',thick);
    set(g, {'color'}, {[0 0.4470 0.7410]; [0.6350 0.0780 0.1840];[0.4660 0.6740 0.1880]});
    hold on
    errorbar(timepoint,TANData{i},TANerr{i},"ok",'LineWidth',thick)
    realizationsA=zeros(length(timepoints),size(chain,1));
    realizationsC=realizationsA;
    realizationsT=realizationsA;
    realizationsTAN=realizationsA;

    variables.I = conditions(i,1);
    variables.T = conditions(i,2);
    variables.u = conditions(i,3);
    c0=1;

    for j = 1:size(chain,1)
        c= chain(j,:);
        paramIN = c;
        realization = plotNumericalInt(paramIN,variables);

        s2 =sqrt(s2chain(j));

        realizationsA(:,j) = realization(:,1)+randn(9,1)*s2;
        realizationsC(:,j)= realization(:,2)+randn(9,1)*s2;
        realizationsT(:,j)= realization(:,3)+randn(9,1)*s2;
        realizationsTAN(:,j)=realizationsA(:,j)+realizationsC(:,j)+realizationsT(:,j);
    end
    
    %Get confidence intervals
    intervalA = prctile(realizationsA,[95,50,5],2);
    intervalC = prctile(realizationsC,[95,50,5],2);
    intervalT = prctile(realizationsT,[95,50,5],2);
    intervalTAN= prctile(realizationsTAN,[95,50,5],2);

    a=shadedErrorBar(timepoint,intervalA(:,2)',[(intervalA(:,1)-intervalA(:,2))';(intervalA(:,2)-intervalA(:,3))'],'lineProps',{'b','LineWidth',thick});
    c=shadedErrorBar(timepoint,intervalC(:,2)',[(intervalC(:,1)-intervalC(:,2))';(intervalC(:,2)-intervalC(:,3))'],'lineProps',{"Color",[0.6350 0.0780 0.1840],'LineWidth',thick});
    t=shadedErrorBar(timepoint,intervalT(:,2)',[(intervalT(:,1)-intervalT(:,2))';(intervalT(:,2)-intervalT(:,3))'],'lineProps',{'Color',[0.4660 0.6740 0.1880],'LineWidth',thick});
    tt=shadedErrorBar(timepoint,intervalTAN(:,2)',[(intervalTAN(:,1)-intervalTAN(:,2))';(intervalTAN(:,2)-intervalTAN(:,3))'],'lineProps',{"k",'LineWidth',thick});
    

    title(["I ="+ conditions(i,1)+ 'mA/cm^2'; "T ="+(conditions(i,2)-273)+ degC])
    xlabel("Time (h)")
    xlim([0,420/60]);
    ylim([0,1.2]);
    ylabel("Non-Dimensional Concentration")
    [h,p] = vartest2([intervalA(:,2),intervalC(:,2),intervalT(:,2)],dataC{i});
    meanFit = [meanFit;h,p];
end
legend("Anode data", "Cathode data", "Trap Data", "Total Data", "Anode Model", "Cathode Model", "Trap Model", "Total Model")
sgtitle("Fitting Set")

%% Validation Set

timepointsSS = [0,30,60,90,120,150]';
variables.timepoints =timepointsSS;

dataset = {ssDNHT,ssDHT,ssINHT, ssIHT};
errset = {errssDNHT,errssDHT,errssINHT, errssIHT};
TANset = {sum(ssDNHT,2),sum(ssDHT,2),sum(ssINHT,2), sum(ssIHT,2)};
TANerrset = {sqrt(sum(errssDNHT.^2,2)),sqrt(sum(errssDHT.^2,2)),sqrt(sum(errssINHT.^2,2)),sqrt(sum(errssIHT.^2,2))};
paramsM = results.mean;
hf=figure;
points=length(timepointsSS);
timepoint=timepointsSS/60;
meanFit = [];
for i = 1:size(conditionsSS,1)
    hf(i)= subplot(2,2,i);
    
    g=errorbar(timepoint,dataset{i},errset{i},"o",'LineWidth',thick);
    set(g, {'color'}, {[0 0.4470 0.7410]; [0.6350 0.0780 0.1840];[0.4660 0.6740 0.1880]});
    hold on
    errorbar(timepoint,TANset{i},TANerrset{i},"ok",'LineWidth',thick)
    variables.I = conditionsSS(i,1);
    variables.T = conditionsSS(i,2);
    variables.u = conditionsSS(i,3);
    realizationsA=zeros(length(timepointsSS),size(chain,1));
    realizationsC=realizationsA;
    realizationsT=realizationsA;
    realizationsTAN=realizationsA;
    for j = 1:size(chain,1)
        c= chain(j,:);
        paramIN = c;
        realization = plotNumericalInt(paramIN,variables);

        s2 =sqrt(s2chain(j));

        realizationsA(:,j) = realization(:,1)+randn(points,1)*s2;
        realizationsC(:,j)= realization(:,2)+randn(points,1)*s2;
        realizationsT(:,j)= realization(:,3)+randn(points,1)*s2;
        realizationsTAN(:,j)=realizationsA(:,j)+realizationsC(:,j)+realizationsT(:,j);
    end
    
    intervalA = prctile(realizationsA,[95,50,5],2);
    intervalC = prctile(realizationsC,[95,50,5],2);
    intervalT = prctile(realizationsT,[95,50,5],2);
    intervalTAN= prctile(realizationsTAN,[95,50,5],2);
   a=shadedErrorBar(timepoint,intervalA(:,2)',[(intervalA(:,1)-intervalA(:,2))';(intervalA(:,2)-intervalA(:,3))'],'lineProps',{'b','LineWidth',thick});
    c=shadedErrorBar(timepoint,intervalC(:,2)',[(intervalC(:,1)-intervalC(:,2))';(intervalC(:,2)-intervalC(:,3))'],'lineProps',{"Color",[0.6350 0.0780 0.1840],'LineWidth',thick});
    t=shadedErrorBar(timepoint,intervalT(:,2)',[(intervalT(:,1)-intervalT(:,2))';(intervalT(:,2)-intervalT(:,3))'],'lineProps',{'Color',[0.4660 0.6740 0.1880],'LineWidth',thick});
    tt=shadedErrorBar(timepoint,intervalTAN(:,2)',[(intervalTAN(:,1)-intervalTAN(:,2))';(intervalTAN(:,2)-intervalTAN(:,3))'],'lineProps',{"k",'LineWidth',thick});
  
    meanFit = [meanFit;[intervalA(:,2),intervalC(:,2),intervalT(:,2)]];
    xlabel("Time (h)")
    ylabel("Non-Dimensional Concentration")
    xlim([0 timepoint(end)])
    ylim([-.2 1.2])
    title(["I ="+ conditionsSS(i,1)+ 'mA/cm^2'; "T ="+(conditionsSS(i,2)-273)+degC])
end
legend("Anode data", "Cathode data", "Trap Data", "Total Data", "Anode Model", "Cathode Model", "Trap Model", "Total Model")
sgtitle("Validation Set")

