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

evolvetime = 10;

nature.evolve(evolvetime);
model.evolve(evolvetime);

naturetomodel = csl.datools.observation.Observation(nature.NumVars);


observeindicies = [2:2:18, 20:40];
nobsvars = numel(observeindicies);
%obserrormodel = csl.datools.error.Tent('NumVars', nobsvars, ...
%    'InitialState', linspace(0.01, 0.49, nobsvars).');

obserrormodel = csl.datools.error.Gaussian('CovarianceSqrt', eye(nobsvars));
observation = csl.datools.observation.Indexed(model.NumVars, ...
    'ErrorModel', obserrormodel, ...
    'Indicies', observeindicies);
R = speye(nobsvars);

ensN = 10;
distfn = model.DistanceFunction;

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