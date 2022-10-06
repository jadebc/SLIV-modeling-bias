% Adaptive MCMC, using Haario et al:
% https://link.springer.com/article/10.1007/s11222-008-9110-y
% and
% http://probability.ca/jeff/ftpdir/adaptex.pdf

% F:          Function giving log-posterior density for a parameter set x
% x0:         Initial value of parameter set x
% n:          Number of iterations
% cov0:       Initial covariance matrix
% fac:        Scaling factor for covariance matrix. Set fac = 1 for default
% fixinds:    Elements of x that should be held fixed. Set fixinds = [] for full MCMC
% blockinds:  Number of 'epi parameters' (e.g. beta, X2) in x, if we want to vary epi and non-epi parameters as independent 'blocks'. Set blockinds = [] if we want a full covariance matrix
% displ:      Structure with display options. Set displ = true to show progress

function [xsto, outsto, history, accept_rate] = MCMC_adaptive(F, x0, n, sigma, fixinds, blockind, cov0, displ)

d = length(x0); b = 0.05; sd = sigma*2.4^2/d;
if ~isempty(fixinds)
    inds = fixinds(1,:); %vals = fixinds(2,:);
    vals = x0(fixinds);
else
    inds = []; vals = [];
end

if isempty(cov0)
    % Checks on the initial covariance matrix
    cov0 = eye(d)*1e-3; cov0(inds,:) = 0; cov0(:,inds) = 0;
    cov0(1:blockind,blockind+1:end) = 0; cov0(blockind+1:end,1:blockind) = 0;
end

% Initiate the output matrices
xsto = zeros(d,n); outsto = zeros(1,n);
history = zeros(d+1,n);                                                    % Rows: 1:d Proposed values 4. Accept or reject

xsto(:,1) = x0(:); xbar = xsto;
%keyboard;
FX = F(x0); outsto(1) = FX;
acc = 0;

if displ; figure; end

% --- Start the MCMC loop -------------------------------------------------
for t = 2:n
    X = xsto(:,t-1);
    
    % --- Make a proposal from the distribution
    Y0 = mvnrnd(X,0.1^2*cov0*sigma/d);
    if t < 2*d
        %Y = max(Y0,0); 
        Y = Y0;
        Y(inds) = vals;
    else
        ind0 = 1; ind1 = t-1;
        covmat = cov(xsto(:,ind0:ind1)');
        covmat(inds,:) = 0; covmat(:,inds) = 0;
        covmat(1:blockind,(blockind+1:end)) = 0; covmat((blockind+1):end,1:blockind) = 0;
        % Need to modify this bit so it goes recursively - faster
        
        % Precaution to make sure values don't get negative
%         if ~all(eigs(sd*covmat)) > eps
%            keyboard; 
%            % mvnrnd(X,sd*covmat)
%         end
        mat = sd*covmat;
        mat = (mat+mat.')/2;
%         keyboard;

        Y = max((1-b)*mvnrnd(X,mat) + b*Y0,0);
        Y(inds) = vals;
    end
    history(1:d,t) = Y;
    
    % --- Decide whether to accept or not
    FY = F(Y);
    if rand < exp(FY-FX)
        % Accept
        xsel = Y(:);
        FX = FY;
        acc = acc+1;
        history(end,t) = 1;
    else
        % Reject
        xsel = xsto(:,t-1);
    end
    xsto(:,t) = xsel;
    outsto(t) = FX;
    xbar(:,t) = (xbar(:,t-1)*(t-1) + xsel)/t;
    
    % Display options
    if displ && (mod(t,round(n/25))==0); fprintf('%0.5g ', t/n*25); end
    if displ && (mod(t,200)==0)
        plot(xsto(setdiff(1:d,fixinds),1:t-1)'); xlim([0 n]); title(sprintf('%0.5g',acc/t));
        drawnow;
    end
    
%     if xsel(5) == 0
%         keyboard;
%     end
end

accept_rate = acc/n;
xsto = xsto';
history = history';