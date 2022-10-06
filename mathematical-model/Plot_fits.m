clear all; load MCMC_res3d;

ix0 = 7e3; nx = 100;
dx = round((size(xsto,1)-ix0)/nx);
xs = xsto(ix0:dx:end,:,1);


% figure; plot(os)
mk = round(size(xs,1)/25); 
for ii = 1:size(xs,1)
    if mod(ii,mk) == 0; fprintf('%0.5g ', ii/mk); end
    [~,aux0] = get_objective(xs(ii,:), p, r, i, s, xi, prm, agg, sel, gps, incid);
    incsto(:,:,ii) = aux0.inc;
end
fprintf('\n');


tis = {'<4yr','5 - 11yo','12 - 17yo','18 - 49yo','50 - 64yo','>65yo'};
inc_prct = permute(prctile(incsto,[2.5,50,97.5],3),[1,3,2]);
figure; fs = 14;

y2 = incid./p.report;
for ii = 1:length(gps.ages)
    subplot(2,3,ii); hold on;
    
    % Plot projections
    y = inc_prct(:,:,ii)';
%     if ii == 4
%        y = y*0.6; 
%     end
%     if ii == 6
%        y = y/0.75; 
%     end
    
    qq = plot(1:length(y), y(2,:));
    jbfill(1:length(y), y(3,:), y(1,:), 'b', 'None', 1, 0.3); hold on;
    
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

