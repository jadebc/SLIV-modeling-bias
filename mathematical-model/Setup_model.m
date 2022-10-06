clear all;

% --- Bring in all the data

load vacc_data;                                                            % Using Get_vaccine_data2
load epid_data;                                                            % Using Get_epi_data
load popn_data;                                                            % Using Get_age_groups3

% Contra Costa selection
opt       = 2;
prm.N     = agenums(:,opt)';
incid_sto = squeeze(data_aligned(:,:,opt,:));                              % Absolute numbers

prm.contact_matrix = new_CMat;


% --- State addresses -----------------------------------------------------

gps.states = {'S', 'Ia','Is','R'};
gps.ages   = {'a1','a2','a3','a4','a5','a6'};
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

r.gamma        = 1/5;
r.beta         = 0.01;
p.sym          = 0.5;                                                      % Prop. of illnesses that are symptomatic
p.relinf_asym  = 0.5;                                                      % Relative infectiousness of asymp
p.report       = p_mode_report;

p.imm  = 0.1*ones(1,length(gps.ages));
p.VE   = VE.mode;
p.vacc = VC.mode;

names = {'r_beta','r_dur','p_sym','p_relinf_asym','p_imm','log_pI0','p_report'};
lgths =        [1,      1,      1,              1,      6,        1,         6];

xi = []; lim = 0;
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    xi.(names{ii}) = inds;
    lim = inds(end);
end
xi.nx = lim;

bds = zeros(xi.nx,2);
bds(xi.r_beta,:)        = [0 0.2];
bds(xi.r_dur,:)         = [1 15];
bds(xi.p_sym,:)         = [0.1 0.9];
bds(xi.p_relinf_asym,:) = [0.1 0.9];
bds(xi.p_imm,:)         = [zeros(length(gps.ages),1), ones(length(gps.ages),1)];
bds(xi.log_pI0,:)       = [1 6];
bds(xi.p_report,:)      = [p.report*0.5; p.report*1.5]';

prm.bounds = bds';
save Model_setup3;
