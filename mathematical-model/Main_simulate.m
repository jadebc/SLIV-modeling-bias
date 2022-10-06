clear all;

% --- Bring in all the data

load vacc_data;
load epid_data;
load popn_data;

% Contra Costa selection
opt   = 2;
% prm.N = (agenums(:,opt)/sum(agenums(:,opt)))';
% incid = data_aligned(:,:,opt,end)/sum(agenums(:,opt))*1e3;               % Rates per 1,000 population
prm.N = agenums(:,opt)';
incid = data_aligned(:,:,opt,end);                                         % Absolute numbers

% figure; incid(20:end,:) = nan; incid = incid./p_mode_report;
% for ii = 1:6
%    subplot(2,3,ii);
%    plot(incid(:,ii)*1e3);
% end

prm.contact_matrix = new_CMat;
% mat = data_aligned(:,:,:,end)./permute(repmat(prm.N,[1,52,2]),[2,1,3]);
% mat = squeeze(sum(data_aligned(:,:,:,end),2))./sum(agenums,1);
% figure; plot(mat)



% --- State addresses -----------------------------------------------------

gps.ages   = {'a1','a2','a3','a4','a5','a6'};
gps.states = {'S', 'Ia','Is','R'};
gps.vaccs  = {'v0','v1'};

[i, s, d, lim] = get_addresses({gps.states, gps.ages, gps.vaccs}, [], [], [], 0);
d = char(d);
i.aux.inc = lim + [1:6];


% --- Counter for age-specific, symptomatic incidence ---------------------

m = zeros(i.nstates);
m(s.Is,:) = 1; 
sel.inc = sparse(m - diag(diag(m)));

m = zeros(length(gps.ages),i.nstates);
for ii = 1:length(gps.ages)
   m(ii,intersect(s.Is, s.(gps.ages{ii}))) = 1;
end
agg.inc = sparse(m);


% --- Parameter specification ---------------------------------------------

r.gamma            = 1/5;
r.beta             = 0.01;

p.sym              = 0.5;                                                  % Prop. of illnesses that are symptomatic
p.relinf_asym      = 0.5;                                                  % Relative infectiousness of asymp

p.report           = p_mode_report;

p.imm  = 0.1*ones(1,length(gps.ages));
p.VE   = VE.mode;
p.vacc = VC.mode;


names = {'r_dur','r_beta','p_sym','p_relinf_asym','p_imm','log_pI0'};
lgths =       [1,       1,      1,              1,      6,        1];
% Also: VE, proportion vaccinated, proportion reported

xi = []; lim = 0;
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    xi.(names{ii}) = inds;
    lim = inds(end);
end
xi.nx = lim;

bds = zeros(xi.nx,2);
bds(xi.r_dur,:)         = [1 15];
bds(xi.r_beta,:)        = [0 0.2];
bds(xi.p_sym,:)         = [0.1 0.9];
bds(xi.p_relinf_asym,:) = [0.1 0.9];
bds(xi.p_imm,:)         = [zeros(length(gps.ages),1), ones(length(gps.ages),1)];
bds(xi.log_pI0,:)       = [1 6];
prm.bounds = bds';






% --- Do random samples and find best starting points ---------------------

nsam = 1e4;
xsam = prm.bounds(1,:) + repmat(diff(prm.bounds,1),nsam,1).*lhsdesign(nsam,xi.nx);

F = @(x) get_objective(x, p, r, i, s, xi, prm, agg, sel, gps, incid);
mk = round(nsam/25);
outs = [];
tic
for ii = 1:nsam
    if mod(ii,mk)==0; fprintf('%0.5g ',ii/mk); end
    xs = xsam(ii,:);
    outs(ii) = F(xs);    
end
toc
mat  = sortrows([outs; 1:nsam]',-1);
ord  = mat(:,2);
ordx = xsam(ord,:);
ordo = outs(ord);
fprintf('\n');









% --- Do fminsearch to get best fit
% x0 = [4.3896    0.0259    0.8315    0.8920    0.0001    0.3777    0.5180    0.8208    0.2459    0.0087    3.4587];
x0 = ordx(1,:);
F = @(x) -get_objective(x, p, r, i, s, xi, prm, agg, sel, gps, incid);
options = optimset('PlotFcns',@optimplotfval);
[x1, fval] = fminsearch(F, x0, options);
[out, aux] = get_objective(x1, p, r, i, s, xi, prm, agg, sel, gps, incid);



% --- Calculate R0
M = aux.M;
% Order of infected indices
inds = [intersect(s.Ia,s.v0), intersect(s.Is,s.v0), intersect(s.Ia,s.v1), intersect(s.Is,s.v1)];

m = zeros(i.nstates);
templ = M.lambda;
v0num = (prm.N.*(1-p.imm).*(1-p.vacc))';
v1num = (prm.N.*(1-p.imm).*(p.vacc))';

% Asymptomatic, incl vaccine effect on susceptibility
m(intersect(s.Ia,s.v0),:) = (1-p.sym)*templ.*v0num;
m(intersect(s.Ia,s.v1),:) = (1-p.sym)*templ.*v1num.*(1-p.VE');
% Symptomatic, incl vaccine effect on susceptibility
m(intersect(s.Is,s.v0),:) = p.sym*templ.*v0num;
m(intersect(s.Is,s.v1),:) = p.sym*templ.*v1num.*(1-p.VE');


% Layers for F matrix
vec  = ones(1,length(gps.ages));                                           % Branching into symptom status
l1   = [(1-p.sym)*vec, p.sym*vec, (1-p.sym)*vec, p.sym*vec];
vec  = (1-p.imm).*prm.N;                                                   % Initial conditions, by prior immunity and vaccine status 
l2   = [(1-p.vacc).*vec, (1-p.vacc).*vec, p.vacc.*vec, p.vacc.*vec];
vec  = ones(1,length(gps.ages));                                           % Reduced susceptibility by vaccination
l3   = [vec, vec, (1-p.VE).*vec, (1-p.VE).*vec];
mat1 = repmat((l1.*l2.*l3)',1,length(inds));
mat2 = repmat(M.lambda(:,inds),4,1);
F    = mat1.*mat2;

V    = -M.lin(inds,inds);
R0   = max(eigs(F*inv(V));






return;

x0 = [5, 0.01, 0.5, 0.5, 0.1*ones(1,6), 4];
% x0 = [11.4659    0.0033    0.8947    0.1917    0.2250    0.1788    0.0000    0.0668    0.1122    0.0008    1.3351];

[out, aux] = get_objective(x0, p, r, i, s, xi, prm, agg, sel, gps, incid);

return;

% F = @(x) get_objective(x, p, r, i, s, xi, prm, agg, sel, gps, incid);
% [xsto, outsto, history, accept_rate] = MCMC_adaptive(F, x0, 1e4, 1, [], [], [], 1);
% fprintf('\n');

inds = find(outsto==max(outsto));
x1 = xsto(inds(1),:);

x1 = ordx(1,:);

[out, aux] = get_objective(x1, p, r, i, s, xi, prm, agg, sel, gps, incid);
y1 = aux.inc;
y2 = incid./p.report;
figure; 
for ii = 1:6
   subplot(2,3,ii); hold on;
   plot(y1(:,ii));
   plot(y2(:,ii));
end

