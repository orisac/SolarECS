function [t,A_t] = extrapolateNumericalInt(params,variables,scenario)
%scenario is a function defining temperature and current conditions at time
%t
%Params come from the DRAM algorithm
%Variables set outside the script
a1 = params(1);
a2 = params(2);
a3 = params(3);
a4 = params(4);
a5 = params(5);

u = variables.u;

c0=1;
timepoints =variables.timepoints;
T_B = @(t) scenario('T',t); %Catholyte Bottle temp
T_C = @(t) 0.42980*T_B(t) + 171.14; %Cathode chamber temp
T_A = @(t) 0.6143*T_C(t) + 114.59; %Anode chamber temp

k1 = @(t) a1*scenario('I',t)*exp(a2/T_A(t));
k2 = @(t) a3*T_C(t)^0.5*exp(a4/T_C(t));
k3 = @(t) a5*T_B(t)^0.5*exp(a4/T_B(t));

dAdt = @(t,A) [-k1(t)*A(1)+scenario('Q',t)*(c0-A(1)),k1(t)*A(1)-(k2(t)+k3(t))*A(2),k2(t)*A(2)]';

[t,A_t] = ode15s(dAdt,timepoints,[c0,0,0]);

end
