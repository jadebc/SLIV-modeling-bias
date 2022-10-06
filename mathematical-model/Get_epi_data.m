clear all;

[num,txt,raw] = xlsread('data.xlsx');

% Organise the incidence data

ctys = {'OUSD','WCCUSD'};
for ic = 1:2
    
    cty = ctys{ic};
    inds     = find(contains(raw(:,4),cty));
    mat      = cell2mat(raw(inds,[2,3,6,7]));
    mat(:,1) = mat(:,1)-2010;
    %mat(:,2) = mod(mat(:,2)-39,52);
    
    nyrs = max(mat(:,1));
    nwks = max(mat(:,2));
    nage = max(mat(:,3));
    
    dat = zeros(nyrs,nwks,nage);
    for iyr = 1:nyrs
        for iwk = 1:nwks
            for ia = 1:nage
                ind = intersect(intersect(find(mat(:,1)==iyr), find(mat(:,2)==iwk)), find(mat(:,3)==ia));
                if ~isempty(ind)
                    dat(iyr,iwk,ia) = mat(ind,4);
                end
            end
        end
    end
    % Combine weeks 52 and 53
    dat(:,end-1,:) = dat(:,end-1,:) + dat(:,end,:);
    dat(:,end,:) = [];
    
    % Record the outputs
    data1(:,:,:,ic) = dat;
    data2(:,:,:,ic) = cat(2,dat(1:end-1,40:end,:), dat(2:end,1:39,:));
end

data_raw     = permute(data1,[2,3,4,1]);                                   % Dims: 1.Week, 2.Age gp, 3.City, 4.Year
data_aligned = permute(data2,[2,3,4,1]);

save epid_data data_raw data_aligned;

% Data_raw is just the existing data
% Data_aligned is the data reorganised, so the epidemic starts at time zero