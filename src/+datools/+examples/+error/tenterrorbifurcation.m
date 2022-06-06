initialevolve = 1000;
addpoints = 2000;
nmu = 2000;

mus = linspace(1, 2, nmu);

xs = zeros(nmu, addpoints);
ys = zeros(nmu, addpoints);

for mui = 1:nmu

    mu = mus(mui);

    tenterror = csl.datools.error.Tent('Mu', mu, 'Scale', 1);

    for i = 1:initialevolve
        tenterror.adderr(0, 0);
    end

    for i = 1:addpoints
        y = tenterror.adderr(0, 0.5);
        xs(mui, i) = mu;
        ys(mui, i) = y;
    end

end

p = scatter(xs(:), ys(:), 1, '.k');
p.MarkerEdgeAlpha = 0.1;
