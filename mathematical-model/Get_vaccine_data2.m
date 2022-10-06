
[num, txt, raw] = xlsread('Vaccine_data.xlsx','Sheet2');

% --- Vaccine coverage numbers (17/18 only) -------------------------------
[r,c] = find(strcmp(raw,'Vaccine coverage'));
VC.mode = cell2mat(raw(r+[2:7],c+2))'/100;

% --- Vaccine efficacy numbers (17/18 only) -------------------------------
[r,c] = find(strcmp(raw,'1718'));
mat = max(cell2mat(raw(r+1:r+5,c-1:c+1)),0);

% Fit beta distribution to these intervals
tmp = [];
for ir = 1:size(mat,1)
    [logfn, tmp(ir,:), aux] = get_distribution_fns(mat(ir,:)/100, 'beta', 1);
end
% Allow an extra row for the 5-11/12-17 split
VE.beta = [tmp(1:2,:); tmp(2:end,:)];

tmp = mat(:,1)'/100;
VE.mode = [tmp(1:2), tmp(2:end)];


% --- Multipliers for each age group --------------------------------------

vec = [143.44, 364.71, 178.16, 94.3, 11];
[num, txt, raw] = xlsread('Multipliers.xlsx');

[r,c] = find(cellfun(@(x)isequal(x,1718),raw));
mat = cell2mat(raw(r,unique(c)+[1:3]))

% Get appropriate log normal distribution parameters
tmp = [];
for ir = 1:size(mat,1)
    [logfn, tmp(ir,:), aux] = get_distribution_fns(mat(ir,:)*vec(ir), 'lognorm', 0);
end
% Allow an extra row for the 5-11/12-17 split
mult_lognorm = [tmp(1:2,:); tmp(2:end,:)];
% Save the peak values
tmp = mat(:,1)'.*vec;
p_mode_report = 1./[tmp(1:2), tmp(2:end)];


save vacc_data VE VC mult_lognorm p_mode_report;


% 0  - 4
% 5  - 11
% 12 - 17
% 18 - 49
% 50 - 64
% 65+