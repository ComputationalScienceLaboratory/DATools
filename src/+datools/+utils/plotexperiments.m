function plotexperiments(filepath)
% load the experiments
load(filepath);

% define the number of figures
f1 = figure;
f2 = figure;

switch filtertype
    case 'Ensemble'
        infs = infs;
    case 'Particle'
        infs = rejs;
end

ensNsplot = ensNs(rankhistogramplotindex);
infsplot = infs(rankhistogramplotindex);
totalrunsplot = numel(ensNsplot) * numel(infsplot);

% find relevant indices
[plotindices, count] = findindices(rankhistogramplotindex,rankhistogramplotindex, numel(ensNs));

if (count ~= totalrunsplot)
    error('Mismatch in plotting index. Check the function findindices \n');
end

for runn = 1:totalrunsplot
    rw = numel(infsplot) - 1 - floor((runn - 1)/numel(ensNsplot));
    cl = runn - floor((runn - 1)/numel(ensNsplot)) * numel(ensNsplot);
    row = floor((runn - 1)/numel(ensNsplot)) + 1;
    col = runn - (row - 1) * numel(ensNsplot);

    ensN = ensNsplot(col);

    ylabelposition = -0.095;
    subxlabelposition = 0.045;
    subylabelposition = 0.1;

    figure(f1);
    ha = subplot(numel(infsplot), numel(ensNsplot), rw*numel(ensNsplot)+cl);
    hold all;
    hold all;
    z = rankvalmatrix{plotindices(runn)};
    maxz = max(z);
    z = z / sum(z);
    NN = numel(z);
    z = NN * z;
    h = bar(xvalmatrix{plotindices(runn)}, z, 'hist');
    h.EdgeColor = 'none';
    h.FaceColor = '#0072BD';
    plot(xvalmatrix{plotindices(runn)}, polyvalmatrix{plotindices(runn)}, 'Color', '#A2142F', 'LineWidth', 3);
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
            ylabelpointer = ylabel(han, 'Inflation', 'FontSize', 10, 'FontWeight', 'bold');
            ylabelpointer.Position = [ylabelposition, 0.5000, -7.1054e-15];
        case 'Particle'
            ylabelpointer = ylabel(han, 'Rejuvetion', 'FontSize', 10, 'FontWeight', 'bold');
            ylabelpointer.Position = [ylabelposition, 0.5000, -7.1054e-15];
    end
    xlabel(han, 'Ensemble Size', 'FontSize', 10, 'FontWeight', 'bold');
    title(han, 'Rank Histogram');
    if (row == 1)
        hapos = get(ha, 'position');
        a = annotation('textbox', [hapos(1) + hapos(3) / 2 - 0.030, subylabelposition, 0, 0], 'string', num2str(ensNsplot(col)));
        a.FontWeight = 'demi';

    end
    if (col == 1)
        hapos = get(ha, 'position');
        a = annotation('textbox', [subxlabelposition, hapos(2) + hapos(4) / 2 + 0.02, 0, 0], 'string', num2str(infsplot(row)));
        a.FontWeight = 'demi';
    end


    figure(f2);
    subplot(numel(infsplot), numel(ensNsplot), rw*numel(ensNsplot)+cl);
    plot(spinup+1:1:steps, rmsvalmatrix{plotindices(runn)});
    xlim([spinup + 1, steps]);
    ylim([0, variance]);
    set(gca, 'XTick', [spinup + 1, steps]);
    set(gca, 'XTickLabel', [spinup + 1, steps]);
    set(gca, 'YTick', [0, variance]);
    set(gca, 'YTickLabel', [0, variance]);
    han = axes(f2, 'visible', 'off');
    han.Title.Visible = 'on';
    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';
    ylabelpointer = ylabel(han, 'Value', 'FontSize', 10, 'FontWeight', 'bold');
    ylabelpointer.Position = [ylabelposition, 0.5000, -7.1054e-15];
    xlabel(han, 'Ensemble Size', 'FontSize', 10, 'FontWeight', 'bold');
    title(han, 'RMSE');
    if (row == 1)

        hapos = get(ha, 'position');
        a = annotation('textbox', [hapos(1) + hapos(3) / 2 - 0.030, subylabelposition, 0, 0], 'string', num2str(ensNsplot(col)));
        a.FontWeight = 'demi';
    end
    if (col == 1)
        hapos = get(ha, 'position');
        a = annotation('textbox', [subxlabelposition, hapos(2) + hapos(4) / 2 + 0.02, 0, 0], 'string', num2str(infsplot(row)));
        a.FontWeight = 'demi';
    end

end

f3 = figure;
f4 = figure;

figure(f3);
switch filtertype
    case 'Ensemble'
        imagesc(ensNs, infs, rmses.');
        caxis([0, variance]);
        colorbar;
        set(gca, 'YDir', 'normal');
        axis square;
        title('RMSE Heatmap');
        colormap('pink');
        xlabel('Ensemble Size', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Inflation', 'FontSize', 10, 'FontWeight', 'bold');
        set(gca, 'XTick', linspace(ensNs(1), ensNs(end), size(ensNs, 2)));
        set(gca, 'XTickLabel', ensNs);
        set(gca, 'YTick', linspace(infs(1), infs(end), size(infs, 2)));
        set(gca, 'YTickLabel', infs);
    case 'Particle'
        imagesc(ensNs, rejs, rmses.');
        caxis([0, variance]);
        colorbar;
        set(gca, 'YDir', 'normal');
        axis square;
        title('RMSE HeatMap');
        colormap('pink');
        xlabel('Ensemble Size', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Rejuvetion', 'FontSize', 10, 'FontWeight', 'bold');
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
        title('D_{KL}');
        colormap(map);
        xlabel('Ensemble Size', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Inflation', 'FontSize', 10, 'FontWeight', 'bold');
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
        title('D_{KL}');
        colormap(map);
        xlabel('Ensemble Size', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Rejuvetion', 'FontSize', 10, 'FontWeight', 'bold');
end

end

function [plotindices, count] = findindices(xindices, yindices, rowsize)
count = 0;
for i = 1:numel(xindices)
    for j = 1:numel(yindices)
       count = count+1;
       plotindices(count) = (xindices(i)-1)*rowsize + yindices(j);
    end
end
end