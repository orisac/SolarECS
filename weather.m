function temp = weather(time)
%contains high and low temperature for Jan 1 or Jul 1 in each location


%Set desired season and locations
season = "winter";
location = "UG";
time = time/60;
if location == "PA"
    %san jose monthly avg
    if season == "summer"
        highlow = [15.6,27.2]; %summer
    else
        highlow = [6,15]; %winter
    end
elseif location == "KY"
    if season == "summer"
        highlow = [17.8,25]; %summer - December
    else
        highlow = [15.6,22.2]; %winter - July
    end
    elseif location == "UG"
    if season == "summer"
        highlow = [17.8,26.7]; %summer - December
    else
        highlow = [16.7,25]; %winter - July
    end

else
    %Oklahoma Weather - montly avg
    if season == "summer"
        highlow = [22,34.4]; %summer
    else
        highlow = [-2.2,8.3]; %winter
    end

end


temp = -(diff(highlow)/2)*cos(2*pi*time/24)+mean(highlow)+273; %[12,32] low and high temp
end