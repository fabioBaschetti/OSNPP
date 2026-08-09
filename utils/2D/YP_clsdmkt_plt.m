function YP_clsdmkt_plt(t_path,x1_path,x2_path,P_min,P_max)

    figure()
    hold on
    h1=plot(t_path,x2_path, 'Color', 'red','LineWidth',0.5);
    hold on
    yline(P_max,'Color','k','LineStyle','--','LineWidth',0.5)
    hold on
    yline(P_min,'Color','k','LineStyle','--','LineWidth',0.5)
    hold on
    h2=plot(t_path,x1_path,'Color','blue','LineStyle','-','LineWidth',0.5);
    grid on
    hold on
    yyaxis left
    ylim([0 ; 1.2])
    set(gca,'YColor','red')
    yyaxis right
    set(gca,'YColor','b')
    ylim([0 ; 1.2])
    yticks([P_min P_max])
    yticklabels({'$P_{\min}$','$P_{\max}$'})
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    legend
    hPatch = findobj(gca,'Type','patch'); 
    for k = 1:length(hPatch)
        hPatch(k).Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    h = findobj(gca,'Type','ConstantLine');
    for k = 1:length(h)
        h(k).Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    legend([h1 h2], ...
           {'$Y_t$','$P_t$'}, ...
           'Interpreter','latex', ...
           'Location','northwest', ...
           'FontSize',12)

end
