clc;
nvar = 40;
nvarnature = 2*nvar;

if nvarnature == nvar
    naturetotruth = speye(nvar);
else
    naturetotruth = spdiags(repmat([1/4, 1/4, 1/2, 1/4, 1/4], nvarnature, 1), [-nvarnature+1, -1, 0, 1, nvarnature-1], nvarnature, nvarnature);
    naturetotruth = naturetotruth(1:2:nvarnature, :);
end

% yes, I know about 17
rng(17);

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
time = spinup*11;

% forecast this much
futurefore = 10;
for i = 1:futurefore
    naturefo.forecast();
    da.forecast();
end

err = zeros(time, 1);

% other good number
%net.numInputs = numel(H) + nvar;
%net.numLayers = 4;

%net = train(net, randn(70, 100), randn(30, 100));


obs = zeros(numel(H), time);
modelgood = zeros(nvar, time);
modelf = zeros(nvar, time);

dstr = '';

for i = 1:time
    truth = naturetotruth*naturefo.CurrentBestGuess;
    
    y = truth(H) + Rsqrt*randn(numel(H), 1);
    
    obs(:, i) = y;
    
    da.analysis(y, R);
    
    modelf(:, i) = mean(da.CurrentBestGuess, 2);

    xa = mean(da.CurrentBestGuess, 2);
    
    modelgood(:, i) = xa;
    
    err(i) = rms(xa - truth);
    
    cumerr = sqrt(mean(err((spinup + 1):i).^2));
    
    if (mod(i, 10) == 0)
    fprintf(repmat('\b', 1, numel(dstr)));
    dstr = sprintf('i: %d, err: %.5f\n', i, cumerr);
    
    fprintf(dstr);
    end
    
    
    % forecast
    naturefo.forecast();
    da.forecast();
        
end

numlayer = 5;

%net = feedforwardnet((numel(H) + nvar)*ones(1, numlayer));
net = feedforwardnet((numel(H))*ones(1, numlayer));



% netinput = [obs(:, (spinup + 2):(end - 1)); modelgood(:, (spinup + 1):(end - 2))];
% netoutput = obs(:, (spinup + 3):end);

%netinput = [obs(:, (spinup + 2):(end - 1)); modelgood(:, (spinup + 2):(end - 1))];
%netoutput = obs(:, (spinup + 3):end);
%netoutput = modelgood(H, (spinup + 3):(end));


netinput = [modelgood(:, (spinup + 2):(end - 1)); obs(:, (spinup + 2):(end - 1)) - modelgood(H, (spinup + 2):(end - 1))];
netoutput = obs(:, (spinup + 3):(end)) - modelgood(H, (spinup + 3):(end));


%netoutput = obs(:, (spinup + 3):end);
%netoutput = modelgood(H, (spinup + 3):(end));

%netoutput = modelgood(:, (spinup + 3):end);

net.trainParam.epochs = 10000;
net.trainParam.max_fail = 1000;

net = train(net, netinput, netoutput, 'useGPU','yes');

save('lorenz96bigmodelobsforenet2.mat', 'net');


return;


errnet = zeros(time, 1);

fprintf('\n\n');
dstr = '';


model.Y0 = da.Problem.Y0;
valfore = csl.datools.DAmethod(model, solver);


valtime = 200;

% for i = 1:valtime
%     truth = naturetotruth*naturefo.CurrentBestGuess;
%     
%     y = truth(H) + Rsqrt*randn(numel(H), 1);
%     
%     %da.analysis(y, R);
%     
%     xa = net([y; valfore.CurrentBestGuess]);
%     
%     valfore.Problem.Y0 = xa;
%     
%     err(i) = rms(xa - truth);
%     
%     cumerr = sqrt(mean(err(1:i).^2));
%     
%     fprintf(repmat('\b', 1, numel(dstr)));
%     dstr = sprintf('i: %d, err: %.5f\n', i, cumerr);
%     fprintf(dstr);
%     
%     
%     % forecast
%     naturefo.forecast();
%     valfore.forecast();
%     
% end
