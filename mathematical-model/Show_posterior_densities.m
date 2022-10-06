clear all;

% load MCMC_1718_out_const;
load MCMC_1718_res_uncon;
xs1 = xsto(:,:,2);

load MCMC_1718_out_const;
xs2 = xsto(:,:,2);

names = {'\beta, symptomatic infection','\gamma, rate of recovery','p_{sym}, proportion symptomatic',{'d, relative infectiousness','symptomatic vs asymptomatic'},{'p_{susc}, proportion','susceptible, 0-4yo'},{'p_{susc}, proportion','susceptible, 5-11yo'},{'p_{susc}, proportion','susceptible, 12-17yo'},{'p_{susc}, proportion','susceptible, 18-49yo'},{'p_{susc}, proportion','susceptible, 50-64yo'},{'p_{susc}, proportion','susceptible, >65yo'},'Initial infected','Proportion reported, 0-4yo','Proportion reported, 5-11yo','Proportion reported, 12-17yo','Proportion reported, 18-49yo','Proportion reported, 50-64yo','Proportion reported, >65yo'};

ff=figure;
for ii = 1:xi.nx
    subplot(6,3,ii); hold on;
    histogram(xs1(:,ii),'Normalization','probability');
    histogram(xs2(:,ii),'Normalization','probability');
    vec = prm.bounds(:,ii);
    %line(xlim,1/diff(vec)*[1 1],'linestyle','--');
    title(names{ii});
end

set(ff,'Position',[440 1 691 796]);


return;

chain = 1;

% Get projections for reported incidence
tmp1    = permute(repmat(xs(:,xi.p_report,chain),[1,1,size(incsto,1)]),[3,2,1]);
mat     = incsto(:,:,:,chain).*tmp1;
inc_plt = permute(prctile(mat,[2.5,50,97.5],3),[3,1,2]);                   % Dims: 1.Lo/Md/Hi, 2.Time, 3.Age

tis = {'<4yr','5 - 11yo','12 - 17yo','18 - 49yo','50 - 64yo','>65yo'};

ff=figure; fs = 14;

% cols = linspecer(size(inc_pct,4));
cols = linspecer(2);

y2 = incid;
for ii = 1:length(gps.ages)
    subplot(2,3,ii); hold on;
    
    % Plot data
    pp = plot(1:20,y2(1:20,ii),'r');
    plot(20:size(y2,1), y2(20:end,ii),'--','Color',pp.Color);
    
    % Plot projections
    y = inc_plt(:,:,ii);
    qq = plot(1:size(y,2), y(2,:), 'b'); hold on;
    jbfill(1:size(y,2), y(3,:), y(1,:), 'b', 'None', 1, 0.3); hold on;
    
    xlim([0 30]);
    set(gca,'fontsize',fs);
    title(tis{ii});
    
    if ismember(ii,[1,4])
        ylabel('Weekly incidence');
    end
    if ismember(ii,4:6)
        xlabel('Week');
    end
    if ii == 2
        legend([pp;qq],'Data','Model projections');
    end
    
end

set(ff,'Position',[680   528   711   450]);

% Get the percent reduction by vaccination
vec = prctile(impact(:,end),[2.5,50,97.5],1)'*100;
fprintf('Reduction in >65 hosps by vaccination: %0.3g (%0.3g - %0.3g) \n\n', vec(2), vec(1), vec(3));


























% inc_pct_uncon = inc_pct;
% save inc_pct_uncon inc_pct_uncon;
