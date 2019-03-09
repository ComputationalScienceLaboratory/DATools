
model = csl.odetestproblems.qgso.presets.GC('large');

[t, y] = ode45(model.F, [0 100], model.Y0);

model.Y0 = y(end, :).';

h = 5;
tend = 5000;
[t, y] = ode45(model.F, 0:h:tend, model.Y0);

VNm1 = y(1:(end - 1), :).';
VN   = y(2:end, :).';

k = 200;
[U, Sigma, W] = svds(VNm1, k);

Stilde = U'*VN*W/Sigma;

dmdprop = @(v) U*(Stilde*(U'*v));

ds = dmdprop(VNm1) - VN;

rms(ds(:))


save('qgdmd.mat', 'dmdprop');

