clear all; 

load Model_setup2;
load Store/MCMC_1718_res_const;      xs1 = xs;
load Store/MCMC_1718_res_uncon;    xs2 = xs;

xrat = xs1./xs2;

mat = [xrat(:,:,1); xrat(:,:,2); xrat(:,:,3); xrat(:,:,4); xrat(:,:,5)];

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