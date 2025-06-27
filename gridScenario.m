function quant = gridScenario(QoI,t,idens,variables)
%Q = .47*1e-3;
Q = variables.Q;

%electrolysis curve, Shen et al
KrE0=[473.674755979412	0.106768024579491	2.30937300460694];
E0 = KrE0(3);
electrolysis = @(V) (V+2*KrE0(1)*KrE0(2).*(V-KrE0(3))-sqrt(V.^2+4*KrE0(1)*KrE0(2)*KrE0(3).*(V-KrE0(3))))/(2*KrE0(2)*(1+KrE0(1)*KrE0(2)));
gridcurrent = idens;
solveV= @(V) electrolysis(V)-gridcurrent;
Vgrid = fzero(solveV ,4);

T = 273+22; %everything at room temp


if QoI == "I"
    quant = gridcurrent;
elseif QoI == "T"
    quant = T;
elseif QoI == "Q"
    quant = Q;
elseif QoI == "P" %power from the grid
    quant =Vgrid*gridcurrent*.064;
end
