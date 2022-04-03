initialevolve = 2000;
addpoints = 2000;
nr = 2000;

rs = linspace(3, 4, nr);

xs = zeros(nr, addpoints);
ys = zeros(nr, addpoints);

for mui = 1:nr

    r = rs(mui);

    tenterror = datools.error.Logistic('R', r, 'Scale', 1);

    % initial evolve
    for i = 1:initialevolve
        tenterror.adderr(0, 0);
    end

    % add points
    for i = 1:addpoints
        y = tenterror.adderr(0, 0.5);
        xs(mui, i) = r;
        ys(mui, i) = y;
    end

end

p = scatter(xs(:), ys(:), 1, '.k');
p.MarkerEdgeAlpha = 0.1;
