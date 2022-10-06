function [p,r,prm] = alloc_parameters(x,p,r,prm,xi)

r.beta             = x(xi.r_beta);
r.gamma            = 1/x(xi.r_dur);
p.sym              = x(xi.p_sym);
p.relinf_asym      = x(xi.p_relinf_asym);
p.imm              = x(xi.p_imm);
p.I0               = 10^(-x(xi.log_pI0));
p.report           = x(xi.p_report); 

% vec = [1, x(xi.cont_matr)];
% prm.contact_matrix = reshape(vec,6,6);