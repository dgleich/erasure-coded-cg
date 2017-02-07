% ufgetmats
% Get the UF matrices for the conditioning study
addpath('~/data/uf');
UFget('refresh');
%%
index = UFget;
mats = find(index.nrows >= 500 & index.ncols <= 10000 & index.posdef & index.numerical_symmetry);
%%
for mat=mats
    UFget(mat);
end
%%
M = length(mats);
kfracs = [0.01; 0.05; 0.1; 0.2];
K = length(kfracs);
symmat = @(A) triu(A,1) + triu(A,1)' + diag(diag(A));
nrep = 25;

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
        parfor t=1:nrep
            E = randn(n,k)*1/sqrt(n);
            A2 = [A A*E; E'*A E'*A*E];
            lam2 = sort(real(eig(full(symmat(A2))))); 
            lam2(abs(lam2) < n*eps(max(lam2))) = []; % remove null-space
            kappa2 = max(lam2)/min(lam2);
            fprintf('n = %8i, kf = %.2f  kappa = %18g kappa2 = %18g, ratio=%.2f\n', n, kfrac, kappa, kappa2, kappa2/kappa);
            results(mi,ki,t) = kappa2/kappa;
        end
    end
end
toc(t0);    

%% If these didn't get computed above
kappas = zeros(M,1);
parfor mi=1:M    
    Prob = UFget(mats(mi));
    A = Prob.A;
    kappas(mi) = cond(full(A));
end

%%
names = index.Name(mats);
groups = index.Group(mats);
save conduf.mat results mats names groups kfracs nrep M kappas

%%  
load conduf
index = UFget;
ns = index.nrows(mats)';
%%
medresults = median(results,3);
plot(ns(1:M),medresults','.');
xlabel('n');
ylabel('\kappa''/\kappa');


%% Filter out experiments with high condition numbers
subset = kappas.*ns < 5e15 & kappas >= 50.;
medresults = median(results(subset,:,:),3);
subns = ns(subset);
subkap = kappas(subset);
submats = mats(subset);
maxresults = max(results(subset,:,:),[],3);
minresults = min(results(subset,:,:),[],3);
subnames = names(subset);
subgroups = groups(subset);

%%
%semilogx(subkap.*subns.*(index.nnz(submats)./index.nrows(submats))',medresults','.')
semilogx(subns,medresults','.')
%%
medresults(:,2)./medresults(:,1)
%%
medresults(:,3)./medresults(:,2)

%%
boxplot([medresults(:,1) maxresults(:,1) medresults(:,4),maxresults(:,4)])
%%
[~,p] = sort(maxresults(:,4));
clf;
%plot([medresults(p,1) maxresults(p,1) medresults(p,4),maxresults(p,4)],'.-','LineWidth',0.5,'MarkerSize',18)
ph = plot([maxresults(p,1) maxresults(p,4)],'.-','LineWidth',0.5,'MarkerSize',18);
hold on;
plot(minresults(p,1),'Color',get(ph(1),'Color')*0.5 + 0.5*[1,1,1])
plot(minresults(p,4),'Color',get(ph(2),'Color')*0.5 + 0.5*[1,1,1])
for i=p'
    h = text(i,maxresults(p(i),4),sprintf(' \\kappa=%8.1e, n=%i, %s',subkap(p(i)), ns(p(i)), 'test'));
    set(h,'Rotation',90,'FontSize',7)
end
ylabel('kappa');
set(gca,'XTick',[]);
hold off;
set_figure_size([6,3]);
print('conduf.eps','-depsc2');

%%
semilogx(subkap,[maxresults(:,1) maxresults(:,4)],'.');
%%
clf;
semilogx(subkap,[maxresults(:,4)],'.');
hold on;
%semilogx(subkap,[maxresults(:,4)],'.','Color',[         0    0.4470    0.7410]);
for i=1:size(maxresults,1)
    line(subkap(i)*[1,1], [minresults(i,4) maxresults(i,4)],'Color',[         0    0.4470    0.7410]);
    if maxresults(i,4) > 1.4
        adjust = 0;
        if ns(i) == 3600
            adjust = -0.01;
        end
        h = text(subkap(i),maxresults(i,4)+adjust,...
            sprintf(' n=%i, %s/%s', ns(i), subgroups{i}, subnames{i}));
        set(h,'Rotation',35,'FontSize',7,'Interpreter','none')
        
    end

    %line(subkap(i)*[1,1], [minresults(i,1) maxresults(i,1)],'Color',[0.9290    0.6940    0.1250]);
end
set(gca,'XScale','log');
xlabel('\kappa');
ylabel('\kappa''/\kappa');
box off;
xlim([1e1,1e12]);
set_figure_size([6,4]);
print('conduf-20p.eps','-depsc2');

%%
clf;
semilogx(subkap,[maxresults(:,1)],'.');
hold on;
%semilogx(subkap,[maxresults(:,4)],'.','Color',[         0    0.4470    0.7410]);
for i=1:size(maxresults,1)
    line(subkap(i)*[1,1], [minresults(i,1) maxresults(i,1)],'Color',[         0    0.4470    0.7410]);
    if maxresults(i,1) > 1.2
        adjust = 0;
        if ns(i) == 3600
            adjust = -0.01;
        end
        h = text(subkap(i),maxresults(i,1)+adjust,...
            sprintf(' n=%i, %s/%s', ns(i), subgroups{i}, subnames{i}));
        set(h,'Rotation',35,'FontSize',7,'Interpreter','none')
        
    end

    %line(subkap(i)*[1,1], [minresults(i,1) maxresults(i,1)],'Color',[0.9290    0.6940    0.1250]);
end
set(gca,'XScale','log');
xlabel('\kappa');
ylabel('\kappa''/\kappa');
box off;
xlim([1e1,1e12]);
set_figure_size([6,3]);
print('conduf-1p.eps','-depsc2');


%%
clf;
semilogx(subkap,[maxresults(:,1)],'.');
hold on;
%semilogx(subkap,[maxresults(:,4)],'.','Color',[         0    0.4470    0.7410]);
for i=1:size(maxresults,1)
    line(subkap(i)*[1,1], [minresults(i,1) maxresults(i,1)],'Color',[         0    0.4470    0.7410]);
    %line(subkap(i)*[1,1], [minresults(i,1) maxresults(i,1)],'Color',[0.9290    0.6940    0.1250]);
end
set(gca,'XScale','log');
%%
plot(diff(medresults'))

%% Print out the data for the table
for i=1:size(maxresults,2)
    fprintf('%i%% ', 100*kfracs(i));
    fprintf(' & %.3f ', prctile(maxresults(:,i),50));
    fprintf(' & %.3f ', prctile(maxresults(:,i),80));
    fprintf(' & %.3f ', prctile(maxresults(:,i),95));
    fprintf(' & %.3f ', max(maxresults(:,i)));
    fprintf('\\\\ \n');
end    
    
%% Print off all the matrix names
for i=1:length(subgroups)
    fprintf('\"%s/%s\", \n', subgroups{i}, subnames{i});
end


