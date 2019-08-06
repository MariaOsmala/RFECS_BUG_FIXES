function [pks locs]=peak_call(prob_dist,peak_dist)
%prob_dist are the prediction scores of a single enhancer region, composed
%of adjacent regions having an enhancer score larger than 0.5
%prob_dist is a vector of enhancer scores between 0.5 and 1
%peak_dist=20
%peak_call(vals1,peak_filt_dist);
[pks locs]=max(prob_dist); %finds the maximum value -> pks
%locs
%locs is the first occurrence of the maximal value

if(length(prob_dist)-locs>peak_dist) %true if locs is far away from the region end length(prob_dist)>locs+peak_dist
    prob_dist2=prob_dist(locs+peak_dist+1:end); %add peak_dist+1 to locs, 1->22, the prob_dist length decreases by 21
    [pks1 locs1]=peak_call(prob_dist2,peak_dist); % find the next maximal value in prob_dist2
    pks=[pks pks1]; %maximal values
    locs=[locs locs+peak_dist+locs1]; %location
    pks1=[];locs1=[];
end

if(locs-1>peak_dist) %true if locs need to be further away from the beginning than peak_dist
    prob_dist2=prob_dist(1:locs-peak_dist-1); %if locs=23, prob_dist(1:2)
    [pks1 locs1]=peak_call(prob_dist2,peak_dist);
    pks=[pks pks1];
    locs=[locs locs1];
    pks1=[];locs1=[];
end


end


%if(test-1>peak_dist)
%   test    
%end
