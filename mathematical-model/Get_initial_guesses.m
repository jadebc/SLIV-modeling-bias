% clear all; load Model_setup3;
% option = {'1718','const'};

season = option{1};
constrained = strcmp(option{2},'const');

if strcmp(season,'1718')
    incid = incid_sto(:,:,end);
    incr_vacc = 0.11;
else
    incid = incid_sto(:,:,end-1);
    incr_vacc = 0.07;
end

F    = @(x) get_objective2(x, p, r, i, s, xi, prm, agg, sel, gps, incid, constrained, incr_vacc);
negF = @(x) -get_objective2(x, p, r, i, s, xi, prm, agg, sel, gps, incid, constrained, incr_vacc);


% --- Do random samples and find best starting points ---------------------

nsam = 1e4;
xsam = prm.bounds(1,:) + repmat(diff(prm.bounds,1),nsam,1).*lhsdesign(nsam,xi.nx);

mk = round(nsam/25);
outs = [];
tic
for ii = 1:nsam
    if mod(ii,mk)==0; fprintf('%0.5g ',ii/mk); end
    xs = xsam(ii,:);
    outs(ii) = F(xs);
end
toc
mat  = sortrows([-outs; 1:nsam]',1);
ord  = mat(:,2);
ordx = xsam(ord,:);
ordo = outs(ord);
fprintf('\n');


% % --- Do fminsearch to get best fit
% % x0 = [4.3896    0.0259    0.8315    0.8920    0.0001    0.3777    0.5180    0.8208    0.2459    0.0087    3.4587];

for ii = 1:10
    fprintf('%0.5g ',ii);
    x0 = ordx(ii,:);
    [x1(ii,:), fval(ii)] = fminsearch(negF, x0, optimset('PlotFcns',@optimplotfval));
    %[out, aux] = get_objective(x1(ii,:), p, r, i, s, xi, prm, agg, sel, gps, incid);
end
fprintf('\n');

mat = sortrows([1:length(fval); fval]',2);
ord = mat(:,1);
x2 = x1(ord,:);

save(['startingpoints_',season,'_',option{2}]);






showfits = 0;
if showfits
    xs = x2([1:5],:); incsto = []; repsto = [];
    for ii = 1:size(xs,1)
        x = xs(ii,:);
        [out, aux] = F(x);
        incsto(:,:,ii) = aux.inc(:,:,1);
        repsto(:,:,ii) = aux.inc(:,:,1).*p.report;
        impact(ii)     = aux.impact;
    end
    fprintf('\n');
    dat = incid;
    
    figure;
    % cols = linspecer(size(incsto,3));
    for ii = 1:size(incsto,2)
        %     col = cols(:,ii);
        mat = squeeze(repsto(:,ii,:));
        
        subplot(2,3,ii); hold on;
        plot(mat);
        plot(dat(:,ii),'r--','linewidth',2);
        xlim([0 52]);
    end
    legend('1','2','3','4','5');
end
