function quant = directgridScenario(QoI,t, tTemp, Tt)
Q = .47*1e-3;
sun = @SolarIrradiance;
%electrolysis curve, Shen et al
KrE0=[473.674755979412	0.106768024579491	2.30937300460694];
E0 = KrE0(3);
electrolysis = @(V) (V+2*KrE0(1)*KrE0(2).*(V-KrE0(3))-sqrt(V.^2+4*KrE0(1)*KrE0(2)*KrE0(3).*(V-KrE0(3))))/(2*KrE0(2)*(1+KrE0(1)*KrE0(2)));
gridcurrent = 15;
solveV= @(V) electrolysis(V)-gridcurrent;
Vgrid = fzero(solveV ,4);

tempfunction = @(ts) interp1(tTemp,Tt,mod(ts,24*60));
Tout = tempfunction(t);
Tpanel = Tout(:,1);
Tcathode = Tout(:,2);


if QoI == "I"
    if sun(t) == 0
        current = 0;
    else
        SolarI = solarPanelOutput(Tpanel,sun(t));
        solveI= @(V) electrolysis(V)*.064-SolarI(V);
        Vintersect = fzero(solveI ,[E0,8]);
        if isnan(Vintersect)
            current = 0;
        else
            current = electrolysis(Vintersect);
        end
    end
    if current < 5
        current = gridcurrent;
    end
    quant = current;
elseif QoI == "T"
    %quant = temp;
    quant = Tcathode;
elseif QoI == "Q"
    quant = Q;
elseif QoI== "ST"
    quant = Tpanel;
elseif QoI == "P" %power from the grid
    if sun(t) == 0
        current =0;
    else
        SolarI = solarPanelOutput(Tpanel,sun(t));
        solveI= @(V) electrolysis(V)*.064-SolarI(V);
        Vintersect = fzero(solveI ,[E0,8]);
        if isnan(Vintersect)
            current = 0;
        else
            current = electrolysis(Vintersect);
            P=0;
        end
    end
    if current < 5
        P = Vgrid*gridcurrent*.064;
    end
    quant = P;
end
end