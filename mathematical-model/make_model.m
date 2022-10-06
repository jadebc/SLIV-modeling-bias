function M = make_model(p, r, i, s, prm, gps)

% Linear transitions matrix
m = zeros(i.nstates);
for iv = 1:length(gps.vaccs)
    for ia = 1:length(gps.ages)
        getid = @(st) i.(st).(gps.ages{ia}).(gps.vaccs{iv});
        Ia = getid('Ia');
        Is = getid('Is');
        R  = getid('R');
        
        % Recovery
        sources = [Ia, Is];
        destin  = R;
        rate    = r.gamma;
        m(destin, sources) = m(destin, sources) + rate;
    end
end
M.lin = sparse(m - diag(sum(m)));


% Nonlinear transitions
for ia = 1:length(gps.ages)
    m = zeros(i.nstates);
    for iv = 1:length(gps.vaccs)
        
        getid = @(st) i.(st).(gps.ages{ia}).(gps.vaccs{iv});
        S  = getid('S');
        Ia = getid('Ia');
        Is = getid('Is');
        
        source  = S;
        destins =      [Ia,    Is];
        rates   = [1-p.sym, p.sym]*(1 - (iv==2)*p.VE(ia));
        m(destins, source) = m(destins, source) + rates';
    end
    M.nlin.(gps.ages{ia}) = sparse(m - diag(sum(m)));
end


% Force-of-infection matrix
m = zeros(length(gps.ages),i.nstates);
for iv = 1:length(gps.vaccs)
    vacc = gps.vaccs{iv};
    
    cols = intersect(s.(vacc),s.Ia);
    m(:,cols) = p.relinf_asym*prm.contact_matrix./repmat(prm.N,length(gps.ages),1);
    
    cols = intersect(s.(vacc),s.Is);
    m(:,cols) = prm.contact_matrix./repmat(prm.N,length(gps.ages),1);
    
end
m = m*r.beta;
M.lambda = sparse(m);
