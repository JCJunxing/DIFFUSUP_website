function [chiresult,kk,Simulation_picked] = chisquare(Length_measured,Measured_profile,error_profile_measured,Length_simulated,Simulation_sorted,include_error)
    % kk is the style number for the chi sqaure 0: sum (o-e)^2/e^2; 1: sum (o-e)^2/var^2;
    % chiresult is the result for chisquare
    % Length_measured is the distance for the spot on measured profile
    % Measured_profile is the measured profile 
    % error_profile_measured is the error for the measured profile 
    % Length_simulated is the distance for the spot on modeling profile
    % Simulation_sorted is the modeled profile 
    % Simulation_picked is the value of the picked spot on the modeling profile. The spot has the same the distance as the measured spot
    % (all composition in one matrix, order the same, set in the code)
    
    Simulation_picked=[];
    
    % find the modeling value which could be picked directly
    [C,ia,ib] = intersect(Length_measured,Length_simulated,'rows');
    for k1 = 1:length(ia)
        Simulation_picked(ia(k1),:) = Simulation_sorted(ib(k1),:);
    end
    
    % Calculating the spots left in picking linear method
    Length_need_calc = setdiff(Length_measured,C);
    for k2 = 1:length(Length_need_calc)
        ia2 = find(Length_measured==Length_need_calc(k2));
        idx = find(Length_simulated<Length_need_calc(k2), 1, 'last');
        Simulation_picked(ia2,:) = (Length_need_calc(k2)-Length_simulated(idx))*(Simulation_sorted(idx+1,:)-Simulation_sorted(idx,:))/(Length_simulated(idx+1)-Length_simulated(idx)) + Simulation_sorted(idx,:);
    end
    
    % Chi-square calculation for the different error settings. (smaller=better)
    if isempty(error_profile_measured) || ~isempty(find(isnan(error_profile_measured) == 1)) || ~isempty(find(~error_profile_measured)) || include_error==0
        kk= 0;
        chiresult = sum((Measured_profile-Simulation_picked).^2./((0.05*Measured_profile).^2),'all')/size(Measured_profile,1)/size(Measured_profile,2);
    else
        kk= 1;
        chiresult = sum((Measured_profile-Simulation_picked).^2./(error_profile_measured.^2),'all')/size(Measured_profile,1)/size(Measured_profile,2);
    end
end
