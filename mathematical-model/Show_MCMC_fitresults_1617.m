clear all; 

load MCMC_out_prev_unconstrained;

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

return;

% Visualise the contact matrix terms
mat = permute(xs(:,xi.cont_matr,:),[1,3,2]);
figure;
for ii = 1:size(mat,3)
    subplot(6,6,ii+1);
    plot(mat(:,:,ii));
end
% For which terms is there least variation across runs?
mat2 = squeeze(mean(mat(end-50:end,:,:),1))';
for ii = 1:size(mat2,1)
   vec = mat2(ii,:);
   vari(ii) = (max(vec)-min(vec))/mean(vec);
end

mat2 = mat(end-50:end,:,:);
for ii = 1:size(mat2,3)
    tmp = mat2(:,:,ii);
    vec = tmp(:); 
    vari(ii) = sqrt(var(vec))/mean(vec);
end

ordmat = sortrows([vari;1:length(vari)]',1)
reshape(0:35,6,6)




return;

for ii = 1:1
    [~,aux0] = get_objective(x2(ii,:), p, r, i, s, xi, prm, agg, sel, gps, incid);
    incsto(:,:,ii) = aux0.inc(:,:,1);
    impact(ii)     = aux0.impact;
end
fprintf('\n');


tis = {'<4yr','5 - 11yo','12 - 17yo','18 - 49yo','50 - 64yo','>65yo'};
inc_plt = permute(incsto,[1,3,2]);
figure; fs = 14;

y2 = incid./p.report;
for ii = 1:length(gps.ages)
    subplot(2,3,ii); hold on;
    
    % Plot projections
    y = inc_plt(:,:,ii);
    qq = plot(1:size(y,1), y);
    %jbfill(1:length(y), y(3,:), y(1,:), 'b', 'None', 1, 0.3); hold on;
    
    % Plot data
    pp = plot(1:20,y2(1:20,ii));
    plot(20:size(y2,1), y2(20:end,ii),'--','Color',pp.Color);
    
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
