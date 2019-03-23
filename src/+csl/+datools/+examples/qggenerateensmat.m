
model = csl.odetestproblems.qgso.presets.GC('large');

[t, y] = ode45(model.F, [0 100], model.Y0);

model.Y0 = y(end, :).';

h = 100;
tend = h*100;
[t, ysamples] = ode45(model.F, 0:h:tend, model.Y0);

ysamples = ysamples(2:end, :).';

save('qgsolargeensemble.mat', 'ysamples');

