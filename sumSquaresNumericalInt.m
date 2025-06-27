function sumSquaresAnalytic = sumSquaresNumericalInt(params, data)
%conditions is rows of [I,T,u] (current, temperature, conductivity)
%unpacking the data
timepoints=data.xdata;
expData=data.udata ;
conditions=data.conditions;


a1 = params(1);
a2 = params(2);
a3 = params(3);
a4 = params(4);
a5 = params(5);

model = [];
n = size(conditions,1);
for i = 1:n
    c0 = 1;
    I = conditions(i,1);
    T = conditions(i,2);
    u = conditions(i,3);

    % Experimentally fitted temperatures
    T_B = T; %Catholyte Bottle temp
    T_C = 0.42980*T_B + 171.14; %Cathode chamber temp
    T_A = 0.6143*T_C + 114.59; %Anode chamber temp

    k1 = a1*I*exp(a2/T_A);
    k2 = a3*T_C^0.5*exp(a4/T_C);
    k3 = a5*T_B^0.5*exp(a4/T_B);

    dAdt = @(t,A) [-k1*A(1),k1*A(1)-(k2+k3)*A(2),k2*A(2)]';
    [t,A_t] = ode15s(dAdt,timepoints,[c0,0,0]); %integrate the equations

    model = [model;A_t]; %store them to compare to experimental data
end

%Catch any integration problems
if any(isnan(model), 'all')
    model = zeros(size(model));
end
if size(model,1)<72
    model = zeros(72,3);
end

diff = model - expData;

%if you want to account for the size of error bars in fitting
if data.errScale ==1
    diff = diff./(data.error/min(min(data.error))); %scale by error bars
end

%Calculated the sum of squares error
diffSq = diff.^2;
sumSquaresAnalytic = sum(diffSq, 'all');
end
