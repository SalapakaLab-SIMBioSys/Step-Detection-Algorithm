% Function to encapsulate the Kerssemakers step detection method & automate
% the manual steps

function [Fit] = call_Kerssemakers(data_raw)
close all
verbosity = 0
% Convert 1 dimensional raw data to have an additional column of counts
data_input(:,2) = data_raw;
data_input(:,1) = 1:1:size(data_raw);

[data, indexes,lijst,properties,initval, Steppedness, selectie] = Steps_Find(data_input, verbosity);

% Determine the appropriate number of steps in the signal 
% Do this by finding the most prominent peak in the chigood/chibad ratio
[pks, locs, w, p] = findpeaks(Steppedness(selectie,2));
[m,n] = max(p);
max_loc = locs(n);
doitforthisstepnumber = Steppedness(selectie(max_loc),3);
%[m,n] = max(Steppedness(selectie,2));
%doitforthisstepnumber = Steppedness(selectie(n),3);

figure()
plot(Steppedness(selectie,2))

[dummy, Steppedness, selectie, Fit] = Steps_Evaluate(data, indexes,lijst, properties,initval,doitforthisstepnumber, verbosity);


end

