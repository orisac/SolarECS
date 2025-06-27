function cost = solarECSEconomics(QoI,input,tvec)
% % Calculates economic values including cost of electricity, electricity
% revenues, and fertilizer revenues


%Set these variables for desired season and location
season = "summer";
location = "PA";

if location == "PA"
    price =546.36; %pacific northwest
elseif location == "KY"
    price = 762.54; %USD/MT urea
elseif location == "UG"
    price = 1.49; %USD/kg urea
else
    price = 540;  %oklahoma
end
fertilizer = price/907185; %price per gram ammonium sulfate https://www.ams.usda.gov/mnreports/ams_3657.pdf

if QoI == "NH3" %give the fertilizer price when input is given in grams of nitrogen/L
    cost = fertilizer.*input.*.5*132.14/(2*14);%convert grams nitrogen to grams ammonium sulfate
    if location == "KY"
        cost = price*1e-6*input*60.06/(2*14);%$/grams urea *urea molar mass/grams of N in urea*grams of N recovered)
    elseif location == "UG"
        cost = price*1e-3*input*60.06/(2*14);
    end
elseif QoI == "W" %input given as a vector of watts
    if season == "summer"
        charge= [0,0.09816, 0.11735,0.14484,0.11735,0.09816,0.09822,0.09822,0.12734,0.09822,0.0569,0.09822,0.09822,0.09822,0.09822,0.12734]*.66; %costs as $/kwh
        if location == "PA"
            %San Jose
            month = 7;
            hour_start =[0,14,16,21,23,24];
        elseif location == "KY"
            hour_start =[0,6,22,24];
            charge=[0,5.7,11.4,5.7]*0.0077;
        elseif location == "UG"
            hour_start =[0,5,18,23,24];
            charge=[0,245.4,387.2,521.3,245.5]/3738.33;
        else
            hour_start = [0,14,19,24];
            charge =[0,0.021062,0.078526,0.021062];
        end
        %tvec is in minutes
        thour = mod(tvec/60,24); %get time as hour of the day
        cost =zeros(3,length(tvec));
        for i = 2:length(thour)
            tnet = (tvec(i)-tvec(i-1))/60; %convert minutes to hours
            for j = 2:length(hour_start)
                if thour(i)<=hour_start(j) %loop until you find the right charge period
                    if thour(i-1)>=hour_start(j-1) %make sure both ends of the interval are in the charge period
                        cost(1,i)=charge(j)*tnet*input(i)*1e-3; %convert to kwH and get cost
                        cost(3,i) = charge(j);
                    else %adjust if there is a portion in the last charge period
                        previouscharge = abs(hour_start(j-1)-thour(i-1)); %portion in last charge time
                        currentcharge = abs(thour(i)-hour_start(j-1)); %portion in current charge time
                        cost(1,i)= input(i)*1e-3*(previouscharge*charge(j-1)+currentcharge*charge(j));
                        cost(3,i) = charge(j);
                    end
                    break
                end
            end
        end
        cost(2,:) = cumsum(cost(1,:)); %get cumulative costs
    else %winter season
        if location == "PA"
            hour_start =[0,16,21,24];
            charge = [0.09822,0.09822,0.12734,0.09822];
        elseif location == "KY"
            hour_start =[0,6,22,24];
            charge=[0,5.7,11.4,5.7]*0.0077;
        elseif location == "UG"
            hour_start =[0,5,18,23,24];
            charge=[0,243.6,385.3,529.4,243.6]/3768.0;
        else
            %Oklahoma data
            hour_start =[0, 24];
            charge =[0,0.019934];
        end

        %tvec is in minutes
        thour = mod(tvec/60,24); %get time as hour of the day
        cost =zeros(3,length(tvec));
        for i = 2:length(thour)
            tnet = (tvec(i)-tvec(i-1))/60; %convert minutes to hours
            for j = 2:length(hour_start)
                if thour(i)<=hour_start(j) %loop until you find the right charge period
                    if thour(i-1)>=hour_start(j-1) %make sure both ends of the interval are in the charge period
                        cost(1,i)=charge(j)*tnet*input(i)*1e-3; %convert to kwH and get cost
                        cost(3,i) = charge(j);
                    else %adjust if there is a portion in the last charge period
                        previouscharge = abs(hour_start(j-1)-thour(i-1)); %portion in last charge time
                        currentcharge = abs(thour(i)-hour_start(j-1)); %portion in current charge time
                        cost(1,i)= input(i)*1e-3*(previouscharge*charge(j-1)+currentcharge*charge(j));
                        cost(3,i) = charge(j);
                    end
                    break
                end
            end
        end
        cost(2,:) = cumsum(cost(1,:)); %get cumulative costs
    end
elseif QoI == "sell" %input given as vector of watts
    if season == "summer"
        %data taken from pge 2023 export prices for jul 1st
        if location == "PA"

            hour_start = 0:24;
            charge =[0,0.03491,0.05138,0.0634,0.0621,0.06338,0.06045,0.0613,0.06125,0.05303,0.04768,0.04779,0.04784,0.04806,0.0486,0.05041,0.05096,0.0522,0.05446,0.0544,0.05336,0.05296,0.05336,0.05294,0.05223];

            %oklahoma data
        elseif location == "KY"
            hour_start =[0,6,22,24];
            charge=[0,5.7,11.4,5.7]*0.0077;
        elseif location == "UG"
            hour_start =[0,5,18,23,24];
            charge=.8*[0,245.4,387.2,521.3,245.5]/3738.33;    
        else
            hour_start = [0,14,19,24];
            charge =[0,0.021062,0.078526,0.021062];
        end
        thour = mod(tvec/60,24); %get time as hour of the day
        cost =zeros(3,length(tvec));
        for i = 2:length(thour)
            tnet = (tvec(i)-tvec(i-1))/60; %convert minutes to hours
            for j = 2:length(hour_start)
                if thour(i)<=hour_start(j) %loop until you find the right charge period
                    if thour(i-1)>=hour_start(j-1) %make sure both ends of the interval are in the charge period
                        cost(1,i)=charge(j)*tnet*input(i)*1e-3; %convert to kwH and get cost
                        cost(3,i) = charge(j);
                    else %adjust if there is a portion in the last charge period
                        previouscharge = abs(hour_start(j-1)-thour(i-1)); %portion in last charge time
                        currentcharge = abs(thour(i)-hour_start(j-1)); %portion in current charge time
                        cost(1,i)= input(i)*1e-3*(previouscharge*charge(j-1)+currentcharge*charge(j));
                        cost(3,i) = charge(j);
                    end
                    break
                end
            end
        end
        cost(2,:) = cumsum(cost(1,:)); %get cumulative costs
    else

        if location == "PA"
            hour_start = 0:24;
            charge =[0,0.05966,0.06533,0.06668,0.06664,0.06679,0.06409,0.05848,0.05448,0.05322,0.05374,0.05403,0.05384,0.05352,0.05403,0.05782,0.06048,0.06089,0.05838,0.05808,0.05645,0.05548,0.0545,0.0545,0.05608];
        elseif location == "KY"
            hour_start =[0,6,22,24];
            charge=[0,5.7,11.4,5.7]*0.0077;
        elseif location == "UG"
            hour_start =[0,5,18,23,24];
            charge=.8*[0,243.6,385.3,529.4,243.6]/3768.0;
        else
            %Oklahoma data
            hour_start =[0, 24];
            charge =[0,0.019934];
        end

        thour = mod(tvec/60,24); %get time as hour of the day
        cost =zeros(3,length(tvec));
        for i = 2:length(thour)
            tnet = (tvec(i)-tvec(i-1))/60; %convert minutes to hours
            for j = 2:length(hour_start)
                if thour(i)<=hour_start(j) %loop until you find the right charge period
                    if thour(i-1)>=hour_start(j-1) %make sure both ends of the interval are in the charge period
                        cost(1,i)=charge(j)*tnet*input(i)*1e-3; %convert to kwH and get cost
                        cost(3,i) = charge(j);
                    else %adjust if there is a portion in the last charge period
                        previouscharge = abs(hour_start(j-1)-thour(i-1)); %portion in last charge time
                        currentcharge = abs(thour(i)-hour_start(j-1)); %portion in current charge time
                        cost(1,i)= input(i)*1e-3*(previouscharge*charge(j-1)+currentcharge*charge(j));
                        cost(3,i) = charge(j);
                    end
                    break
                end
            end
        end
        cost(2,:) = cumsum(cost(1,:)); %get cumulative costs
    end
end

end
