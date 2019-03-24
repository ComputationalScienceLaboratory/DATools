
model = otp.qg.presets.PopovSandu2019('Nature');

[t, y] = ode45(model.Rhs.F, [0 100], model.Y0);

model.Y0 = y(end, :).';

h = 100;
tend = h*250;
[t, ysamples] = ode45(model.Rhs.F, 0:h:tend, model.Y0);

ysamples = ysamples(2:end, :).';

save('qglargeensemble.mat', 'ysamples');

