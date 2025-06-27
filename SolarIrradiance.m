function irradiance = SolarIrradiance(time)
%Loads Solar Profiles

%set desired curve
season = "winter";
location = "UG"

time = time/60;
% %summer
if location == "PA"
    if season == "summer"
PaloAltoJul1sun=readmatrix("PaloAltoJul1sun.csv");
sun=@(t) interp1([0;PaloAltoJul1sun(:,1);24],[0;PaloAltoJul1sun(:,2);0]*1e3,mod(t,24));
irradiance=sun(time);
% 
% %winter
    else
PaloAltoJan1sun=readmatrix("PaloAltoJan1sun.csv");
sun=@(t) interp1([0;PaloAltoJan1sun(:,1);24],[0;PaloAltoJan1sun(:,2);0]*1e3,mod(t,24));
irradiance=sun(time);
    end

%Uganda weather
elseif location == "UG"
    if season == "summer"
%summer
PaloAltoJul1sun=readmatrix("UGJan1sun.csv");
sun=@(t) interp1([0;PaloAltoJul1sun(:,1);24],[0;PaloAltoJul1sun(:,2);0]*1e3,mod(t,24));
irradiance=sun(time);
% 
    else
%winter
PaloAltoJan1sun=readmatrix("UGJul1sun.csv");
sun=@(t) interp1([0;PaloAltoJan1sun(:,1);24],[0;PaloAltoJan1sun(:,2);0]*1e3,mod(t,24));
irradiance=sun(time);
    end

    %OK Weather
else
    if season == "summer"
%summer
PaloAltoJul1sun=readmatrix("OKJul1sun.csv");
sun=@(t) interp1([0;PaloAltoJul1sun(:,1);24],[0;PaloAltoJul1sun(:,2);0]*1e3,mod(t,24));
irradiance=sun(time);
% 
    else
%winter
PaloAltoJan1sun=readmatrix("OKJan1sun.csv");
sun=@(t) interp1([0;PaloAltoJan1sun(:,1);24],[0;PaloAltoJan1sun(:,2);0]*1e3,mod(t,24));
irradiance=sun(time);
    end
end
end