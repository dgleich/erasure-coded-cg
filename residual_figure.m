function residual_figure(matname,resvec,failiter,nfail)

clf;
semilogy(resvec,'.-','LineWidth',1);
if failiter > 0
    hold on;
    plot(failiter,resvec(failiter),'ro','MarkerSize',10);
    plot(failiter,resvec(failiter),'rx','MarkerSize',10);
    hold off;
end
ylabel('residual norm');
xlabel('iteration');
xlim([1,numel(resvec)]);
set_figure_size([2.5,2.5])
set(gca,'LineWidth',0.6);
box off;
print(sprintf('%s-%i.eps',matname,nfail),'-depsc2');