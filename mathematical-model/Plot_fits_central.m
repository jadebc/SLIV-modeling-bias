% Plot fits using central estimates gotten from Get_initial_guesses

clear all; load res_sto5;

% figure; plot(os)
% for ii = 1:size(x2,1)
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

