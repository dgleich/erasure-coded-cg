ns = 250:250:3000;
N = length(ns);
%N = 3;
nrep = 25;
kfracs = [0.01; 0.05; 0.1; 0.2];
K = length(kfracs);
symmat = @(A) triu(A,1) + triu(A,1)' + diag(diag(A));

results = zeros(N,K,nrep);

t0 = tic;
for ni=1:N
    n = ns(ni);
    A = gallery('tridiag',n,-1,2,-1);
    kappa = cond(full(A));
    for ki=1:K
        kfrac = kfracs(ki);
        k = ceil(kfrac*n);
        for t=1:nrep
            %Q = qr(randn(n,k),0);
            E = randn(n,k)*1/sqrt(n);
            A2 = [A A*E; E'*A E'*A*E];
            lam2 = eig(full(symmat(A2))); 
            lam2(lam2 < n*eps(max(lam2))) = []; % remove null-space
            kappa2 = max(lam2)/min(lam2);
            fprintf('n = %8i, kf = %.2f  kappa = %18g kappa2 = %18g, ratio=%.2f\n', n, kfrac, kappa, kappa2, kappa2/kappa);
            results(ni,ki,t) = kappa2/kappa;
        end
    end
end
toc(t0);

%%
save condscaling.mat results N ns kfracs

%%
load condscaling
medresults = median(results,3);
maxresults = max(results,[],3);
minresults = min(results,[],3);
hs = plot(ns(1:N),medresults','.-','LineWidth',1.5,'MarkerSize',15);
colors = get(hs,'Color');
clf;
hold on;
for ki=1:length(kfracs)
    color = colors{ki};
    color = color*0.5 + 0.5*[1,1,1];
    p = fill([ns(1:N),fliplr(ns(1:N))],[maxresults(:,ki); flipud(minresults(:,ki))]', color);
    set(p,'EdgeColor','none');
    hs(ki) = plot(ns(1:N),medresults(:,ki),'.-','LineWidth',1.5,'MarkerSize',15,'Color',colors{ki});
    text(ns(end),medresults(end,ki),sprintf(' %i%%',100*kfracs(ki)));
end
hold off;
xlabel('n');
ylabel('\kappa''/\kappa');
xlim([0,3250]);
set_figure_size([3,3]);
box off;

%legend(hs,arrayfun(@(x) sprintf('%i%%',x),100*kfracs,'UniformOutput',false),'Orientation','Horizontal')
%legend boxoff;

print(gcf,'scaling-cond.eps','-depsc2','-painters');