function switches_plt(i0,T,ttau_,iota_,m)

    figure
    stairs([0.00 ttau_ T],[i0 iota_ iota_(end)], 'LineWidth',1.5);
    grid on
    xlabel('t')
    x_shift = 0.025*T;
    xlim([-x_shift,T+x_shift])
    yticks(1:m);
    ylabels = arrayfun(@(j) sprintf('regime : %d', j), 1:m, 'UniformOutput', false);
    yticklabels(ylabels);
    y_shift = 0.025*(m-1);
    ylim([1-y_shift,m+y_shift])
    title(['initial regime: i(0) = ', num2str(i0)])

end