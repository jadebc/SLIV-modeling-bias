% clear all; load startingpoints_prev_unconstrained;

co0 = cov(x2);

tic
for ii = 1:3
    F = @(x) get_objective(x, p, r, i, s, xi, prm, agg, sel, gps, incid, constrained, incr_vacc);
    [outx, outo, history, accept_rate] = MCMC_adaptive(F, x2(ii,:), 2e4, 1, [], [], co0, 1);
    xsto(:,:,ii) = outx;
    outsto(:,:,ii) = outo;
    fprintf('\n');
end
toc

save(['MCMC_',season,'_res_',option{2}]);
beep; pause(0.5); beep; pause(0.5); beep;   