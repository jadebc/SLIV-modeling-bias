function R0 = find_R0(p, r, i, s, prm, gps)

M = make_model(p, r, i, s, prm, gps);

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
% Assumble into F
F = m(inds,inds);

V = -M.lin(inds,inds);
R0 = max(eigs(F*inv(V)));
