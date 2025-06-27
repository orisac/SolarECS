function A_t = plotNumericalInt(params,variables)

a1 = params(1);
a2 = params(2);
a3 = params(3);
a4 = params(4);
a5 = params(5);

T = variables.T;
I = variables.I;

c0=1;
timepoints =variables.timepoints;
T_B = T; %Catholyte Bottle temp
T_C = 0.42980*T_B + 171.14; %Cathode chamber temp
T_A = 0.6143*T_C + 114.59; %Anode chamber temp

k1 = a1*I*exp(a2/T_A);
k2 = a3*T_C^0.5*exp(a4/T_C);
k3 = a5*T_B^0.5*exp(a4/T_B);

dAdt = @(t,A) [-k1*A(1),k1*A(1)-(k2+k3)*A(2),k2*A(2)]';

[t,A_t] = ode15s(dAdt,timepoints,[c0,0,0]);

%Catch any integration problems S
if any(isnan(A_t), 'all') 
    A_t = zeros(length(timepoints),3);
end 
if size(A_t,1)<length(timepoints)
    A_t = zeros(length(timepoints),3);
end
end