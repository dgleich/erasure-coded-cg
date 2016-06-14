%% This experiment tests the number of iterations involved
% as we scale a given matrix. 

ns = 250:250:3000;
N = length(ns);
nrep = 25;
kfracs = [0.01; 0.05; 0.1; 0.2];
K = length(kfracs);
symmat = @(A) triu(A,1) + triu(A,1)' + diag(diag(A));

niters0 = zeros(N,1);
results = zeros(N,K,nrep);

tol=1e-10;

t0 = tic;
for ni=1:N
    n = ns(ni);
    A = gallery('tridiag',n,-1,2,-1);
    xtrue = randn(n,1);
    b = A*xtrue;
    maxit = 10*n;
    x0 = [];
    fail_ind = [];          % start with no failure.
    fail_point = 0;
    num_fail = 0;
    
    [x2_sol_fail, flag, iter, fail_ind, x_fail_ind, resvec] = cg_erasure(...
        A, b, tol, maxit, x0, fail_ind, 1, fail_point, num_fail);
    niters0(ni) = iter;
    
    for ki=1:K
        kfrac = kfracs(ki);
        k = ceil(kfrac*n);        
        for t=1:nrep
            E = randn(n,k)*1/sqrt(n);
            A2 = [A A*E; E'*A E'*A*E];
            b2 = [b; E'*b];
            
            [x2_sol_fail, flag, iter2, fail_ind, x_fail_ind, resvec] = cg_erasure(...
                A2, b2, tol, maxit, x0, fail_ind, 1, fail_point, k);
            results(ni,ki,t) = iter2/iter;
        end
    end
end
toc(t0);

%%
save numiter_scaling.mat results N ns kfracs niters0

%%
load numiter_scaling
medresults = median(results,3);
maxresults = max(results,[],3);
minresults = min(results,[],3);
hs = plot(ns(1:N),medresults','.-','LineWidth',1.5,'MarkerSize',15);
colors = get(hs,'Color');
clf;

ests = [1.45 1.6 1.7 1.85]; % estimates from the scaling plot

hold on;
for ki=1:length(kfracs)
    color = colors{ki};
    color = color*0.5 + 0.5*[1,1,1];
    p = fill([ns(1:N),fliplr(ns(1:N))],[maxresults(:,ki); flipud(minresults(:,ki))]', color);
    set(p,'EdgeColor','none');
    hs(ki) = plot(ns(1:N),medresults(:,ki),'.-','LineWidth',1.5,'MarkerSize',15,'Color',colors{ki});
    text(ns(end),medresults(end,ki),sprintf(' %i%%',100*kfracs(ki)));
    %line(xlim, sqrt(ests(ki))*[1,1]);
end
hold off;
xlabel('n');
ylabel('Iterations-on-Encoded / Iterations-on-original');
xlim([0,3250]);
set_figure_size([3,3]);
box off;




%legend(hs,arrayfun(@(x) sprintf('%i%%',x),100*kfracs,'UniformOutput',false),'Orientation','Horizontal')
%legend boxoff;

print(gcf,'scaling-numiter.eps','-depsc2','-painters');