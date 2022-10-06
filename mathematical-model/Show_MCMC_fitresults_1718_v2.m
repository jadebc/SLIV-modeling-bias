% Version finally used for manuscript: plot in terms of reported cases (to
% stay in same units as data), and show additional outputs for cumulative
% reported cases

clear all;

% load MCMC_1718_out_const;
load MCMC_1718_out_uncon;

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
        ylabel('Weekly reported cases');
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


cinc     = squeeze(sum(mat,1));
cinc_pct = prctile(cinc,[2.5,50,97.5],2)';
cincid   = sum(incid,1);

xx = 1:length(cincid);
figure; hold on; fs = 14; lw = 1.5; ms = 30;
md   = cinc_pct(2,:);
hilo = diff(cinc_pct,1);

cols = {'b','r'};
indslist = {1:3,4:6};

for ii = 1:2
    inds = indslist{ii};
    if ii == 2
        simcol = 'b'; datcol = 'r';
    else
        simcol = 0.5*ones(1,3); datcol = simcol;
    end
    lg(2,:) = plot(xx(inds), md(inds), '.','markersize',ms,'Color',simcol);
    errorbar(xx(inds), md(inds), hilo(1,inds), hilo(2,inds),'linestyle','None','linewidth',lw,'Color',simcol);
    lg(1,:) = plot(xx(inds)+0.1,cincid(inds),'d','markersize',8,'linewidth',lw,'Color',datcol);
    xlim([0.5 6.5]);
end
set(gca,'fontsize',fs,'XTicklabel',tis);
xlabel('Age group');
ylabel('Cumulative reported cases');
legend(lg,'Data','Scaled model projections','Location','NorthWest');












% inc_pct_uncon = inc_pct;
% save inc_pct_uncon inc_pct_uncon;
