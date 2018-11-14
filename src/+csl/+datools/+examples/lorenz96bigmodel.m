clc;
nvar = 40;
nvarnature = 2*nvar;

if nvarnature == nvar
    naturetotruth = speye(nvar);
else
    naturetotruth = spdiags(repmat([1/4, 1/4, 1/2, 1/4, 1/4], nvarnature, 1), [-nvarnature+1, -1, 0, 1, nvarnature-1], nvarnature, nvarnature);
    naturetotruth = naturetotruth(1:2:nvarnature, :);
end



nature = csl.odetestproblems.lorenz96.presets.RandomIC(nvarnature);
model = csl.odetestproblems.lorenz96.presets.RandomIC(nvar);
model.TimeSpan = [0, 0.05];
nature.TimeSpan = [0, 0.05];


ntol = 1e-8;
nJvp = nature.JacobianVectorProduct;

solver      = @(f,t,y) csl.utils.rk4(f, t, y, 0.05);
solvernature = @(f, t, y) csl.utils.RODAS4(f, t, y, ntol, nJvp);

distfn = model.DistanceFunction;

radius = 1.5;
inflation = 1.01;


H = [2:2:(nvar/2 - 1), nvar/2:nvar];

E = speye(nvar);

memain = 1;
meoff  = 1/2;

Qest = spdiags(repmat([meoff, meoff, memain, meoff, meoff], nvar, 1), [-nvar+1, -1, 0, 1, nvar-1], nvar, nvar);

ensN = 10;

parallel = false;
ripits = [1 1];

R = speye(numel(H));
Rsqrt = sqrtm(R);


locfun = csl.datools.statistical.ensemble.localisation.gauss_tiny(radius, distfn);

naturefo = csl.datools.DAmethod(nature, solvernature);

da = csl.datools.statistical.ensemble.DEnKF(model, solver, H,ensN, E, Qest, locfun, inflation, parallel, ripits);

spinup = 500;
time = 5500;

% forecast this much
futurefore = 10;
for i = 1:futurefore
    naturefo.forecast();
    da.forecast();
end

err = zeros(time, 1);

dstr = '';

for i = 1:time
    truth = naturetotruth*naturefo.CurrentBestGuess;
    
    y = truth(H) + Rsqrt*randn(numel(H), 1);
    
    da.analysis(y, R);
    
    err(i) = rms(mean(da.CurrentBestGuess, 2) - truth);
    
    cumerr = sqrt(mean(err((spinup + 1):i).^2));
    
    fprintf(repmat('\b', 1, numel(dstr)));
    dstr = sprintf('i: %d, err: %.5f\n', i, cumerr);
    fprintf(dstr);
    
    
    % forecast
    naturefo.forecast();
    da.forecast();
    
end