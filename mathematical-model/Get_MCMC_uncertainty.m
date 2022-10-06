% Taking the results from the 'free contact matrix' model - simulate results 
% for epidemic dynamics and proportion protection, ready for plotting

clear all; load MCMC_res3d;

ix0 = 1e4; nx = 250; dx = round((size(xsto,1)-x0)/nx);
xs = xsto(ix0:dx:end,:,:);

incsto = []; impact = []; mk = round(size(xs,1)/25);

for iz = 1:3
for ii = 1:size(xs,1)
    if mod(ii,mk)==0; fprintf('%0.5g ',ii/mk); end
    [~,aux0] = get_objective(xs(ii,:,iz), p, r, i, s, xi, prm, agg, sel, gps, incid);
    incsto(:,:,ii,iz) = aux0.inc(:,:,1);
    impact(ii,iz)     = aux0.impact;
end
fprintf('\n');
end

imp_pct = prctile(impact,[2.5,50,97.5],1);
inc_pct = permute(prctile(incsto,[2.5,50,97.5],3),[3,1,2,4]);

save MCMC_res3d_outputs;