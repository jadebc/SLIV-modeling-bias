function out = goveqs(t, in, M, i, agg, sel)

invec = in(1:i.nstates);
out   = zeros(length(invec)+6,1);

lam   = M.lambda*invec;
M_all = M.lin + lam(1)*M.nlin.a1 + lam(2)*M.nlin.a2 + lam(3)*M.nlin.a3 + lam(4)*M.nlin.a4 + lam(5)*M.nlin.a5 + lam(6)*M.nlin.a6;

% Linear terms
out(1:i.nstates) = M_all*invec;
out(i.aux.inc)   = agg.inc*((sel.inc.*M_all)*invec);
