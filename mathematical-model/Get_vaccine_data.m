

% --- Get the vaccine programme numbers
[num, txt, raw] = xlsread('Vaccine_data.xlsx','Vaccine_Data');

VE_dat = cell2mat(raw(4:9,9:11));
VC_dat = cell2mat(raw(14:19,9:10));

% --- Multipliers


% --- Get contact matrix data
[num, txt, raw] = xlsread('Vaccine_data.xlsx','Contact_matrix');
CMat = cell2mat(raw(4:18,3:17));
aggr = cell2mat(raw(23:28,2:16));
new_CMat = aggr*CMat*aggr';


% find(strcmp(raw,'Aggregations'))

% 0  - 4
% 5  - 11
% 12 - 17
% 18 - 49
% 50 - 64
% 65+

return;

% Vaccine coverage estimates - from https://www.cdc.gov/flu/fluvaxview/reportshtml/reporti1718/reportii/index.html
vc_pt = [67.8, 58.3, 48.7, 25.3, 39.3, 58.1];
vc_ui = [6.1, 4.9, 5.4, 2, 2.9, 3.3];

% Vaccine efficacy measures - from https://www.cdc.gov/flu/vaccines-work/2017-2018.html



% Converting data with arbitrary age boundaries, to the boundaries we want
agelims =   [0 4  12  17   49   64  100];
agenums = [809 950 596 4380 2263 1840];
data    = [67.8 58.3 48.7 25.3 39.3 58.1; 6.1 4.9 5.4 2.0 2.9 3.3];

tgt_agelims = [0 4 11 17 49 64 100];

% out = interp1(agelims, data(1,:), tgt_agelims)

for ii = 1:(length(tgt_agelims)-1)
   % Iterate over each window
   lo = tgt_agelims(ii);
   hi = tgt_agelims(ii+1);
   
   % Find the windows to translate from
   ind1 = find(agelims<=lo,1,'last');
   ind2 = find(agelims<=hi,1,'last');
   
   if ismember(ind2,[ind1, ind1+1])
       % New set is a subset of reference ones - just take data directly
       out(ii) = data(ii);
   end
       
   [ind1, ind2]
    
end