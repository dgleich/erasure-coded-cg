% ufgetmats
% Get the UF matrices for the conditioning study
addpath('~/data/uf');
UFget('refresh');
%%
index = UFget;
mats = find(index.nrows >= 500 & index.ncols <= 5000 & index.posdef & index.numerical_symmetry);
%%
for mat=mats
    UFget(mat);
end
%%
M = length(mats);
kfracs = [0.01; 0.05; 0.1; 0.2];
K = length(kfracs);
symmat = @(A) triu(A,1) + triu(A,1)' + diag(diag(A));
nrep = 1;

results = zeros(M,K,nrep);
ns = zeros(M,1);
kappas = zeros(M,1);

t0 = tic;
for mi=1:M    
    Prob = UFget(mats(mi));
    A = Prob.A;
    n = size(A,1);
    ns(mi) = n;
    kappa = cond(full(A));
    kappas(mi) = kappa;
    for ki=1:K
        kfrac = kfracs(ki);
        k = ceil(kfrac*n);
        for t=1:nrep
            E = randn(n,k)*1/sqrt(n);
            A2 = [A A*E; E'*A E'*A*E];
            lam2 = real(eig(full(symmat(A2)))); 
            lam2(lam2 < n*eps(max(lam2))) = []; % remove null-space
            kappa2 = max(lam2)/min(lam2);
            fprintf('m = %2i, n = %8i, kf = %.2f  kappa = %18g kappa2 = %18g, ratio=%.2f\n', ...
                mi, n, kfrac, kappa, kappa2, kappa2/kappa);
            results(mi,ki,t) = kappa2/kappa;
        end
    end
end
toc(t0);    

%%
medresults = median(results,3);
plot(ns(1:M),medresults','.');
xlabel('n');
ylabel('\kappa''/\kappa');

%%
plot(index.nnz(mats)./index.nrows(mats), medresults','.');

%%
semilogx(kappas(1:M), medresults','.');
%%
loglog(kappas(1:M), medresults'*diag(kappas),'.');
%%
loglog(kappas(1:M), medresults'*diag(kappas)-ones(size(medresults'))*diag(kappas+ns),'.');
ylim([0,1e5]);
%%
subset = medresults(kappas < 1e11,:);
subsetkap = kappas(kappas < 1e11);
subsetns = ns(kappas < 1e11);
%%
Z = (diag(subsetkap)*subset - diag(subsetkap)*ones(size(subset)))./(diag(subsetns)*ones(size(subset)));


%%
hist(
