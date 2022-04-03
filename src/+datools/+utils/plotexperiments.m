function plotexperiments(filepath)
% load the experiments
load(filepath);

% define the number of figures
f1 = figure;
f2 = figure;

for runn = 1:totalruns
    rw = numel(infs) - 1 - floor((runn - 1)/numel(ensNs));
    cl = runn - floor((runn - 1)/numel(ensNs)) * numel(ensNs);
    row = floor((runn - 1)/numel(ensNs)) + 1;
    col = runn - (row - 1) * numel(ensNs);

    figure(f1);
    ha = subplot(numel(infs), numel(ensNs), rw*numel(ensNs)+cl);
    hold all;
    hold all;
    z = rankvalmatrix{runn};
    maxz = max(z);
    z = z / sum(z);
    NN = numel(z);
    z = NN * z;
    bar(xvalmatrix{runn}, z);
    plot(xvalmatrix{runn}, polyvalmatrix{runn}, '-*r');
    set(gca, 'XTick', [xs(1), xs(end)]);
    set(gca, 'XTickLabel', [1, ensN + 1]);
    set(gca, 'YTick', []);
    set(gca, 'YTickLabel', []);
    han = axes(f1, 'visible', 'off');
    han.Title.Visible = 'on';
    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';
    switch filtertype
        case 'Ensemble'
            ylabel(han, 'Inflation');
        case 'Particle'
            ylabel(han, 'Rejuvetion');
    end
    xlabel(han, 'Ensemble Size');
    title(han, 'Rank Histogram');
    if (row == 1)

        hapos = get(ha, 'position');
        a = annotation('textbox', [hapos(1) + hapos(3) / 2 - 0.020, 0.1, 0, 0], 'string', num2str(ensNs(col)));
        a.FontWeight = 'demi';

    end
    if (col == 1)
        hapos = get(ha, 'position');
        a = annotation('textbox', [0.065, hapos(2) + hapos(4) / 2 + 0.02, 0, 0], 'string', num2str(infs(row)));
        a.FontWeight = 'demi';
    end


    figure(f2);
    subplot(numel(infs), numel(ensNs), rw*numel(ensNs)+cl);
    plot(spinup+1:1:steps, rmsvalmatrix{runn});
    xlim([spinup + 1, steps]);
    ylim([0, 1]);
    set(gca, 'XTick', [spinup + 1, steps]);
    set(gca, 'XTickLabel', [spinup + 1, steps]);
    set(gca, 'YTick', [0, 1]);
    set(gca, 'YTickLabel', [0, 1]);
    han = axes(f2, 'visible', 'off');
    han.Title.Visible = 'on';
    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';
    ylabel(han, 'Value');
    xlabel(han, 'Time Step');
    title(han, 'RMSE');
    if (row == 1)

        hapos = get(ha, 'position');
        a = annotation('textbox', [hapos(1) + hapos(3) / 2 - 0.020, 0.1, 0, 0], 'string', num2str(ensNs(col)));
        a.FontWeight = 'demi';
    end
    if (col == 1)
        hapos = get(ha, 'position');
        a = annotation('textbox', [0.065, hapos(2) + hapos(4) / 2 + 0.02, 0, 0], 'string', num2str(infs(row)));
        a.FontWeight = 'demi';
    end

end

f3 = figure;
f4 = figure;

figure(f3);
switch filtertype
    case 'Ensemble'
        imagesc(ensNs, infs, rmses.');
        caxis([0, 1]);
        colorbar;
        set(gca, 'YDir', 'normal');
        axis square;
        title('Rmse HeatMap');
        colormap('pink');
        xlabel('Ensemble Size');
        ylabel('Inflation')
        set(gca, 'XTick', linspace(ensNs(1), ensNs(end), size(ensNs, 2)));
        set(gca, 'XTickLabel', ensNs);
        set(gca, 'YTick', linspace(infs(1), infs(end), size(infs, 2)));
        set(gca, 'YTickLabel', infs);
    case 'Particle'
        imagesc(ensNs, rejs, rmses.');
        caxis([0, 1]);
        colorbar;
        set(gca, 'YDir', 'normal');
        axis square;
        title('Rmse HeatMap');
        colormap('pink');
        xlabel('Ensemble Size');
        ylabel('Rejuvetion');
        set(gca, 'XTick', linspace(ensNs(1), ensNs(end), size(ensNs, 2)));
        set(gca, 'XTickLabel', ensNs);
        set(gca, 'YTick', linspace(rejs(1), rejs(end), size(rejs, 2)));
        set(gca, 'YTickLabel', rejs);
end


bn = bone;
pk = flipud(pink);
figure(f4);
map1 = bn;
map1 = map1(51:2:end-1, :);
map2 = pk;
map = [map1; map2(2:2:end-50, :)];
switch filtertype
    case 'Ensemble'
        imagesc(ensNs, infs, rhplotval.');
        caxis([-0.1, 0.1]);
        colorbar;
        set(gca, 'YDir', 'normal');
        set(gca, 'XTick', linspace(ensNs(1), ensNs(end), size(ensNs, 2)));
        set(gca, 'XTickLabel', ensNs);
        set(gca, 'YTick', linspace(infs(1), infs(end), size(infs, 2)));
        set(gca, 'YTickLabel', infs);
        axis square;
        title('KLDiv');
        colormap(map);
        xlabel('Ensemble Size');
        ylabel('Inflation');
    case 'Particle'
        imagesc(ensNs, rejs, rhplotval.');
        caxis([-0.1, 0.1]);
        colorbar;
        set(gca, 'YDir', 'normal');
        set(gca, 'XTick', linspace(ensNs(1), ensNs(end), size(ensNs, 2)));
        set(gca, 'XTickLabel', ensNs);
        set(gca, 'YTick', linspace(rejs(1), rejs(end), size(rejs, 2)));
        set(gca, 'YTickLabel', rejs);
        axis square;
        title('KLDiv');
        colormap(map);
        xlabel('Ensemble Size');
        ylabel('Rejuvetion');
end

end