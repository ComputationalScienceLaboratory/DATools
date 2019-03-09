natureode = csl.odetestproblems.qgso.presets.GC;
modelode  = csl.odetestproblems.qgso.presets.GC;

natureode.TimeSpan = [0, 5];
modelode.TimeSpan = [0 5];

h = 1;
solvernature = @(f, t, y) csl.utils.rk4(f, t, y, round(diff(t)/h));
solvermodel  = solvernature;

% nature will have no sythetic error
naturesyntherror = csl.datools.error.Error;

nature = csl.datools.Model(natureode, solvernature, ...
    'SynthError', naturesyntherror);
model  = csl.datools.Model(modelode,  solvermodel);

naturetomodel = csl.datools.observation.Observation(nature.NumVars);

observeindicies = round(linspace(1, model.NumVars, 300));
nobsvars = numel(observeindicies);
%obserrormodel = csl.datools.error.Tent('NumVars', nobsvars, ...
%    'InitialState', linspace(0.01, 0.49, nobsvars).');

obserrormodel = csl.datools.error.Gaussian('CovarianceSqrt', speye(nobsvars));
observation = csl.datools.observation.Indexed(model.NumVars, ...
    'ErrorModel', obserrormodel, ...
    'Indicies', observeindicies);
R = speye(nobsvars);

ensN = 25;
distfn = model.DistanceFunction;

% We make the assumption that there is no model error
modelerror = csl.datools.error.Error;

load('qgsolargeensemble.mat');
ensembleGenerator = @(x) ysamples;