initialevolve = 100;
addpoints = 500;
nrho = 3000;
tspandt = 50;

rhos = linspace(0, 200, nrho);

rhoplot = zeros(nrho, addpoints);
xs = zeros(nrho, addpoints);
ys = zeros(nrho, addpoints);
zs = zeros(nrho, addpoints);

for rhoi = 1:nrho
    
    rho = rhos(rhoi);
    
    model = otp.lorenz63.presets.Canonical;
    model.Parameters.rho = rho;
    
    model.TimeSpan = [0, initialevolve] + model.TimeSpan(1);
    
    [~, yi] = ode45(model.Rhs.F, [0 initialevolve], model.Y0);
    
    model.TimeSpan = model.TimeSpan + tspandt;
    model.Y0 = yi(end, :).';
    
    rhoplot(rhoi, :) = rho;
    [~, y] = ode45(model.Rhs.F, ...
        linspace(model.TimeSpan(1), model.TimeSpan(end), addpoints), ...
        model.Y0);
    
    xs(rhoi, :) = y(:, 1).';
    ys(rhoi, :) = y(:, 2).';
    zs(rhoi, :) = y(:, 3).';
    
end

xmin = min(xs(:));
xmax = max(zs(:));
zmin = min(zs(:));
zmax = max(zs(:));

xs = (xs(:) + xmin)/diff([xmax, xmin]);
zs = (zs(:) + xmin)/diff([zmax, zmin]);

ds = sqrt(xs.^2 + zs.^2);

[~, si] = sort(ds);

rhoplot = rhoplot(si);
ys = ys(si);

cmap = jet(length(ys));
p = scatter(rhoplot(:), ys(:), 1, cmap);
p.MarkerEdgeAlpha = 0.25;

cax = gca;

set(cax, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca,'XColor','none')
set(gca,'YColor','none')
set(gcf,'color','white')