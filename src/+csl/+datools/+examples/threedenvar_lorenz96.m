clear; clc;

evolvetime = 10;

natureode = csl.odetestproblems.lorenz96.presets.Canonical;
modelode  = csl.odetestproblems.lorenz96.presets.Canonical;

h = 0.05;
solvernature = @(f, t, y) csl.utils.rk4(f, t, y, round(diff(t)/h));
solvermodel  = solvernature;

% nature will have no sythetic error
naturesyntherror = csl.datools.error.Error;

nature = csl.datools.Model(natureode, solvernature, ...
    'SynthError', naturesyntherror);
model  = csl.datools.Model(modelode,  solvermodel);

nature.evolve(evolvetime);
model.evolve(evolvetime);

naturetomodel = csl.datools.observation.Observation(nature.NumVars);


observeindicies = [2:2:18, 20:40];
nobsvars = numel(observeindicies);
obserrormodel = csl.datools.error.Tent('NumVars', nobsvars, ...
    'InitialState', linspace(0.01, 0.49, nobsvars).');
observation = csl.datools.observation.Indexed(model.NumVars, ...
    'ErrorModel', obserrormodel, ...
    'Indicies', observeindicies);
R = speye(nobsvars);

ensN = 10;
inflation = sqrt(1.1);
radius = 4.0;
distfn = model.DistanceFunction;
localization = @(t, y, H) csl.datools.tapering.gc(t, y, radius, distfn, H);

% We make the assumption that there is no model error
modelerror = csl.datools.error.Error;

ensemblePseudo = csl.datools.error.Tent('InitialState', linspace(0.01, 0.49, model.NumVars).', 'Scale', 10);
ensemble = zeros(model.NumVars, ensN);
 
for j = 1:ensN
    for k = 1:10
        ensemblePseudo.adderr(0, 0);
    end
   ensemble(:, j) = ensemblePseudo.adderr(0, 0);
end
ensembleGenerator = @(x) ensemble;

enkf = csl.datools.hybrid.envar.ThreeDEnVar(model, ...
    'Observation', observation, ...
    'NumEnsemble', ensN, ...
    'ModelError', modelerror, ...
    'EnsembleGenerator', ensembleGenerator, ...
    'Inflation', inflation, ...
    'Localization', localization);

spinup = 100;
times = 11*spinup;

mses = zeros(times - spinup, 1);

rmse = nan;

ps = '';

for i = 1:times
    % forecast
    nature.evolve();
    enkf.forecast();
    
    % observe
    xt = naturetomodel.observeWithoutError(nature.TimeSpan(1), nature.State);
    y = enkf.Observation.observeWithError(model.TimeSpan(1), xt);
    
    % analysis
    enkf.analysis(R, y);
    
    xa = enkf.BestEstimate;
    
    if i > spinup
        mses(i - spinup) = mean((xa - xt).^2);
        rmse = sqrt(mean(mses(1:(i - spinup))));
    end
    
    for kk = 1:numel(ps)
        fprintf('\b');
    end
    
    ps = sprintf('step: %d, rmse: %.5f\n', i, rmse);
    
    fprintf(ps);
end
