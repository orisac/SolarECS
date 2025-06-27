function [t,Tt]=TransferHeat(tspan) %t input is in minutes
%TransferHeat calculates the temperature of the solar panel, catholyte
%bottle for a given time. Q and Tamb are function handles of time


tspan = tspan*60; %convert time to seconds
A = (221-2*3.75)*(257-2*9.44)*1e-6; %[m^2] Area of panel
%Qin = @(t) Q*A; %Convert from W/m^2 to W
Rcell = .33; %[W/k] -- fitted from solar sim data. Resistance between catholyte and ambient
Rpanel_water=.16; %[W/k] %resistance between panel and catholyte
Rpanel_air= .81; %[W/k] %resistance between panel and air
mw = 1.5; %mass of water kg (for all liquid in the cell)
cw = 4184; %J/kg
mp =.357; %mass of solar panel, kg -- taken from website
cp = 700;%J/kg from pavlovic et al
Tamb = @weather;
Q=@SolarIrradiance;
dTdt = @(t,T) [1/(mp*cp), 0;0, 1/(mw*cw)]*([Q(t/60)*A+1/Rpanel_air*Tamb(t/60);1/Rcell*Tamb(t/60)]+[-1/Rpanel_air-1/Rpanel_water, 1/Rpanel_water; 1/Rpanel_water, -1/Rpanel_water-1/Rcell]*T);

[t,Tt] = ode23(dTdt,tspan,[Tamb(tspan(1));Tamb(tspan(1))]);
t=t/60;%convert time back to minutes
