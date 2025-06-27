function quant = mpptScenario(QoI,t,tTemp, Tt)
Q = .47*1e-3;
sun=@SolarIrradiance;
current = 500/64;% same as my circuit

tempfunction = @(ts) interp1(tTemp,Tt,mod(ts,24*60));
Tout = tempfunction(t);
Tpanel = Tout(:,1);
Tcathode = Tout(:,2);


if sun(t) == 0
    mpp = 0;
else
SolarI = solarPanelOutput(Tpanel,sun(t));
mpp = max((0:.01:8).*SolarI(0:.01:8)); %get power output at max power point
end


if QoI == "I"
    quant = current;
elseif QoI == "T"
    quant = Tcathode;
elseif QoI == "Q"
    if current ==0
        quant = 0;
    else
        quant = Q;
    end
elseif QoI== "ST"
    quant = Tpanel;
elseif QoI =="P"
    if mpp <0
    mpp =0;
    end
    quant = mpp;
end
end