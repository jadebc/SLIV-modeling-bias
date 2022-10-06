% v2: Version using in-sheet search terms to more accurately locate relevant data 
% v3: Using new (US-based) contact matrix from Prem et al, extracted from
% Excel file 'MUestimates_all_locations_2' in folder 'contact_matrices_152_countries'

[num,txt,raw] = xlsread('StF_age_sex_counts.xlsx','2017ACS');

% --- Population sizes (note can also use the first sheet for data without
% interpolation)
vec = num(:,1);
mat = zeros(6,length(vec));
mat(1,1) = 1;
mat(2,[2,3])   = [1,2/5];
mat(3,[3,4])   = [3/5,3/5];
mat(4,[4:8])   = [2/5,1,1,1,5/10];
mat(5,[8:10])  = [5/10,1,1,];
mat(6,[11:13]) = 1;
agenums = mat*num(:,[1,2]);

% --- Contact matrix data
[num, txt, raw] = xlsread('contact_matrix_data_USA.xlsx','Data');

[r,c] = find(strcmp(raw,'Contact matrix'));
CMat = cell2mat(raw(r+[3:18],c+[2:17]));

[r,c] = find(strcmp(raw,'Aggregations'));
aggr_mat = cell2mat(raw(r+[2:7],c+[1:16]));

[r,c] = find(strcmp(raw,'Demography'));
tmp = cell2mat(raw(r+1,c+[1:16]));
pop_weight = tmp/sum(tmp);

% Aggregate down the rows
tmpa = aggr_mat*CMat;

% Aggregate across the columns
mat = [];
for ii = 1:size(tmpa,1)
    cols = find(aggr_mat(ii,:));                                           % Columns to be aggregated
    tmp  = pop_weight(cols).*aggr_mat(ii,cols); 
    wts  = (tmp/sum(tmp));                              % Aggregation weights
    mat(:,ii) = sum(tmpa(:,cols).*repmat(wts,size(tmpa,1),1),2);
end
new_CMat = mat;    
    
save popn_data agenums new_CMat