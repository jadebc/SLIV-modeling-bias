clear all; 

load MCMC_1718_out_const;
mat1 = squeeze(sum(incsto,1));

load MCMC_1718_out_uncon;
mat2 = squeeze(sum(incsto,1));

% Get the ratio
rat  = mat1./mat2;
rat_pct = permute(prctile(rat,[2.5,50,97.5],2),[2,1,3]);                   % Dims: 1.Lo/Md/Hi, 2.Age, 3.Chain no.



return;


% load MCMC_1718_out_uncon;

% load inc_pct_uncon;

inc_pct = permute(prctile(incsto,[2.5,50,97.5],3),[3,1,2,4]);              % Dims: 1.Lo/Md/Hi, 2.Time, 3.Age, 4.Chain no.

tis = {'<4yr','5 - 11yo','12 - 17yo','18 - 49yo','50 - 64yo','>65yo'};
inc_plt = inc_pct(:,:,:,1);

ff=figure; fs = 14; 

cols = linspecer(size(inc_pct,4));
cols = linspecer(2);

y2 = incid./p.report;
for ii = 1:length(gps.ages)
    subplot(2,3,ii); hold on;
    
    % Plot data
    pp = plot(1:20,y2(1:20,ii),'r');
    plot(20:size(y2,1), y2(20:end,ii),'--','Color',pp.Color);

    % Plot projections
    for ip = 1:1%size(inc_pct,4)
        y = inc_pct(:,:,ii,ip);
        qq = plot(1:size(y,2), y(2,:), 'b'); hold on;
        jbfill(1:size(y,2), y(3,:), y(1,:), 'b', 'None', 1, 0.3); hold on;
        
%         y = inc_pct_uncon(:,:,ii,ip);
%         qq = plot(1:size(y,2), y(2,:), 'g'); hold on;
%         jbfill(1:size(y,2), y(3,:), y(1,:), 'g', 'None', 1, 0.3); hold on;
    end

    xlim([0 30]);
    set(gca,'fontsize',fs);
    title(tis{ii});
    
    if ismember(ii,[1,4])
        ylabel('Symptomatic incidence');
    end
    if ismember(ii,4:6)
        xlabel('Week');
    end
    if ii == 5
        legend([pp;qq],'Scaled data','Model projections');
    end
    
end

set(ff,'Position',[680   528   711   450]);

prctile(impact(:,end),[2.5,50,97.5],1)'*100

% inc_pct_uncon = inc_pct;
% save inc_pct_uncon inc_pct_uncon;
