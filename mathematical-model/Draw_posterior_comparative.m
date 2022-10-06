clear all; 

lbls = {'Infectivity',{'Avg duration','infec'},'Prop symptom','Rel Infec asympto','Prior immune, 0-4','Prior immune, 5-11','Prior immune, 12-17','Prior immune, 18-49','Prior immune, 50-64','Prior immune, 65+','Prior infected','Prop report, 0-4','Prop report, 5-11','Prop report, 12-17','Prop report, 18-49','Prop report, 50-64','Prop report, 65+'};

lbls = {'Infectivity',{'Average duration','of infectiousness'},'Proportion symptomatic',{'Relative infectiousness,','asymptomatic cases'},...
    'Prior immunity, 0-4','Prior immunity, 5-11','Prior immunity, 12-17','Prior immunity, 18-49','Prior immunity, 50-64','Prior immunity, 65+',...
    'Initially infected', {'Proportion symptomatic','cases reported, 0-4'}, {'Proportion symptomatic','cases reported, 5-11'}, ...
    {'Proportion symptomatic','cases reported, 12-27'}, {'Proportion symptomatic','cases reported, 18-49'}, {'Proportion symptomatic','cases reported, 50-64'}, ...
    {'Proportion symptomatic','cases reported, 65+'}};
    
load Model_setup2;
load Store/MCMC_1718_res_const;    xs1 = xs;
load Store/MCMC_1718_res_uncon;    xs2 = xs;

xsc = [xs1(:,:,1); xs1(:,:,2); xs1(:,:,3); xs1(:,:,4); xs1(:,:,5)];
xsu = [xs2(:,:,1); xs2(:,:,2); xs2(:,:,3); xs2(:,:,4); xs2(:,:,5)];

xsc = [xs1(:,:,1); xs1(:,:,2)];
xsu = [xs2(:,:,2); xs2(:,:,3)];

% ff=figure; 
% subplot(1,2,1); hold on;
% plot(xsc); 
% for ii = 1:size(xs1,3)
%     line(151*ii*[1 1], ylim);
% end
%     
% subplot(1,2,2); hold on;
% plot(xsu);
% for ii = 1:size(xs2,3)
%     line(151*ii*[1 1], ylim);
% end


ff=figure;
for ii = 1:size(xsu,2)
    subplot(6,3,ii); hold on;
    histogram(xsu(:,ii));
    histogram(xsc(:,ii));
    title(lbls{ii});
end
set(ff,'Position',[680   178   579   800]);

return;

% return;
% xrat = xs1./xs2;
% mat = [xrat(:,:,1); xrat(:,:,2); xrat(:,:,3); xrat(:,:,4); xrat(:,:,5)];

mat = xsc./xsu;

% xrat_pct = prctile(xrat,[5,50,95],1);
xrat_pct = prctile(mat,[5,50,95],1);

% --- Hacks to fix Infs and zeros
inds = find(xrat_pct(3,:)==Inf);
xrat_pct(3,inds) = 2*xrat_pct(2,inds)-xrat_pct(1,inds);

inds = find(xrat_pct(1,:)==0);
vec = xrat_pct(1,:); vec(inds) = [];
xrat_pct(1,inds) = min(vec);

% --- Plot the results
ys = log10(xrat_pct);
hilo = diff(ys,1);

% --- Plot the results horizontally
ff = figure; fs = 12;
plot(fliplr(ys(2,:)),1:size(ys,2),'.','markersize',24); hold on;
errorbar(fliplr(ys(2,:)), 1:size(ys,2), fliplr(hilo(1,:)), fliplr(hilo(2,:)), 'horizontal', 'linestyle', 'None', 'linewidth', 1.5);
line([0 0], ylim, 'linestyle', '--');
lbls = fliplr({'Infectivity','Avg duration infec','Prop symptom','Rel Infec asympto','Init immune, 0-4','Init immune, 5-11','Init immune, 12-17','Init immune, 18-49','Init immune, 50-64','Init immune, 65+','Init infected','p report, 0-4','p report, 5-11','p report, 12-17','p report, 18-49','p report, 50-64','p report, 65+'});
set(gca,'fontsize',fs,'YTick',[1:size(ys,2)],'YTickLabel',lbls);
xlabel('log_{10} (Constrained/Unconstrained)');