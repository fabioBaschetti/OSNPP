function YP_openmkt_plt(t_path,x1_path,x2_path,P_min_IT,P_max_IT,P_max_EU,ref_path)

    figure() 
    h1=plot(t_path,x2_path, 'Color', 'red');
    hold on
    yline(P_max_IT,'Color','k','LineStyle','--','LineWidth',0.5)
    hold on
    yline(P_min_IT,'Color','k','LineStyle','--','LineWidth',0.5)
    hold on
    h2=plot(t_path,x1_path,'Color','blue','LineStyle','-','LineWidth',1.5);
    grid on
    hold on
    masks = { ref_path < 0, ...
              ref_path >= 0 & ref_path <= P_max_EU, ...
              ref_path > P_max_EU
            };
    colori = [ 0.75 0.95 0.75;
               1.00 1.00 0.70;
               1.00 0.70 0.70
             ];
    yl = [0.0 1.2];
    for j = 1:numel(masks)
        mask = masks{j};
        for k = 1:length(t_path)-1
            if mask(k)
                h = patch([t_path(k) t_path(k+1) t_path(k+1) t_path(k)], ...
                          [yl(1) yl(1) yl(2) yl(2)], ...
                          colori(j,:), 'EdgeColor','none','FaceAlpha',0.4);
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end
        end
    end
    yyaxis left
    set(gca,'YColor','red')
    yyaxis right
    set(gca,'YColor','b')
    ylim([yl(1) yl(2)])
    yticks([P_min_IT P_max_IT])
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
