%% Show the number of iterations for CG based on the condition number change
% This experiment helps empirically justify that if we change
% the condition number by $\sigma$, then we get $\sqrt{\sigma}$ more
% iterations. And this holds when the condition number $\kappa$ is large
% and $\sigma$ is small.
kappas = logspace(2,10);
tol = 1e-8;
iters = log(tol)./(log((sqrt(kappas)-1)./(sqrt(kappas)+1)));
loglog(kappas,iters)
hold all;
sigma = 1.5;
iters2 = log(tol)./(log((sqrt(sigma*kappas)-1)./(sqrt(sigma*kappas)+1)));
loglog(kappas,iters2)
hold off;
%%
plot((iters2./iters).^2)

%%
l = 100;
clf;
semilogy(((sqrt(kappas)-1)./(sqrt(kappas))).^l)
hold all;
semilogy(((sqrt(sigma*kappas)-1)./(sqrt(sigma*kappas))).^(sqrt(sigma)*l))
hold off;

%%
[log((sqrt(kappas)-1)./(sqrt(kappas)+1))' sqrt(sigma)*log((sqrt(sigma*kappas)-1)./(sqrt(sigma*kappas)+1))']

%%
[log((sqrt(kappas)-1)./(sqrt(kappas)+1))' -2./(sqrt(kappas'))]