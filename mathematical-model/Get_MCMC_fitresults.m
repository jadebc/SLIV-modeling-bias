% clear all; source = 'MCMC_1718_res_uncon.mat';
% load(source);

ix0 = 1e4;
nx  = 150;
dx  = round((size(xsto,1)-ix0)/nx);
xs  = xsto(ix0:dx:end,:,:);

incsto = [];
impact = [];
mk = round(size(xs,3)*size(xs,1)/25); ct = 1;
for iz = 1:size(xs,3)
    for ii = 1:size(xs,1)
        if mod(ct,mk)==0; fprintf('%0.5g ', ct/mk); end
        xsam = xs(ii,:,iz); 
        [~,aux0] = get_objective(xsam, p, r, i, s, xi, prm, agg, sel, gps, incid);
        incsto(:,:,ii,iz) = aux0.inc(:,:,1);
        impact(ii,iz)     = aux0.impact;
        
        [p,r,prm] = alloc_parameters(xsam,p,r,prm,xi);
        R0(ii,iz) = find_R0(p, r, i, s, prm, gps);
        
        ct = ct+1;
    end
end
fprintf('\n');

% save(strrep(source,'res','out'));
save(['MCMC_',season,'_out_',option{2}]);