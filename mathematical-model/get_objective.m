function [out, aux] = get_objective(x, p, r, i, s, xi, prm, agg, sel, gps, incid, constrained, incr_vaccine)

% Checks on boundaries
tmp = [prm.bounds; x]; tmp = diff(tmp([1,3,2],:),1);

if min(tmp(:))<0
    out = -Inf; aux = NaN;

else
    [p1,r1,prm1] = alloc_parameters(x,p,r,prm,xi);
        
    pvacc = repmat(p1.vacc,2,1);
    pvacc(2,2) = pvacc(2,2) + incr_vaccine;
    options = odeset('AbsTol',1e-10,'RelTol',1e-10);

    for ii = 1:size(pvacc,1)
        p1.vacc = pvacc(ii,:);
        
        % --- Construct the initial conditions ----------------------------
        init   = zeros(1,i.nstates+6);
        for iv = 1:length(gps.vaccs)
            if iv == 1
                props = 1-p1.vacc;
            elseif iv == 2
                props = p1.vacc;
            end
            inds = intersect(s.S, s.(gps.vaccs{iv}));
            init(inds) = prm.N.*(1-p1.imm).*props;
            inds = intersect(s.R, s.(gps.vaccs{iv}));
            init(inds) = prm.N.*p1.imm.*props;
        end
        
        % --- Construct the model
        M1 = make_model(p1, r1, i, s, prm1, gps);
        
        % --- Seed initial infection amongst unvaccinated -----------------
        init(i.S.a2.v0)  = init(i.S.a2.v0)*(1-p1.I0);
        init(i.Is.a2.v0) = init(i.Is.a2.v0) + init(i.S.a2.v0)*p1.I0;
        
        init(i.S.a2.v1)  = init(i.S.a2.v1)*(1-p1.I0);
        init(i.Is.a2.v1) = init(i.Is.a2.v1) + init(i.S.a2.v1)*p1.I0;
        
        % --- Simulate the dynamics ---------------------------------------
        [t,soln] = ode15s(@(t,in) goveqs(t, in, M1, i, agg, sel), [0 364], init, options);
        
        tmp = interp1(t,soln,0:7:t(end));
        inc = diff(tmp(:,i.aux.inc),1);
        inc(inc<=0) = 1e-15;
        
        incsto(:,:,ii) = inc;
        tmp            = sum(inc,1).*x(xi.p_report);
        cnhosp(ii)     = tmp(end);
    end
    
    
    % --- Matching against the epi data
    mu  = incsto(:,:,1).*p.report;
    sig = sqrt(incsto(:,:,1).*p.report.*(1-p.report));
    dat = incid;
    mat = -log(sig) - (dat-mu).^2./(2*sig.^2);
    
    tmp = mat(9:20,3:6); weights = [1, 10, 1, 10]; weights = [1 1 1 1];
    tmp2 = tmp.*weights;
%     tmp = mat(1:20,:);
    tmp = tmp2(:); out1 = sum(tmp(~isnan(tmp)));
    
    % --- Matching against the hospitalisations
    pha = 1 - cnhosp(2)/cnhosp(1);
    
    a = 37; b = 99;
    out2 = (a-1)*log(pha) + (b-1)*log(1-pha) - gammaln(a) - gammaln(b) + gammaln(a+b);

    % --- Bring them together
    if constrained
        out = out1 + 10*out2;
    else
        out = out1;
    end
%     keyboard;
    
    % --- Additional outputs
    aux.inc  = incsto;
    aux.impact = pha;
    aux.out = [out1 out2];
    %aux.soln = [soln,t];
    %aux.M    = M;
%     keyboard;
end