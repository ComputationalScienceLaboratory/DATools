%% Create Nature

natureode = csl.odetestproblems.qgso.presets.GC('huge', 'low');
%natureode = csl.odetestproblems.qgso.presets.GC('large', 'high');


natureode.TimeSpan = [0, 5];

solvernature = @(f, t, y) ode45(f, t, y);

% nature will have no sythetic error
naturesyntherror = csl.datools.error.Error;
nature = csl.datools.Model(natureode, solvernature, ...
    'SynthError', naturesyntherror);

% create the nature-to-model operator
nf = sqrt(natureode.NumVars);
nc = (nf - 1)/2;
rc = 1:nc;
if2cs = [rc, rc, rc];
jf2cs = [(2*rc)-1, 2*rc, (2*rc)+1];
vf2cs = [1/4*ones(1, nc), 1/2*ones(1, nc), 1/4*ones(1, nc)];
If2c1D = sparse(if2cs, jf2cs, vf2cs, nc, nf);
If2c = kron(speye(nc), If2c1D) * kron(If2c1D, speye(nf));

naturetomodel = csl.datools.observation.Linear(nature.NumVars, 'H', If2c);

%naturetomodel = csl.datools.observation.Observation(nature.NumVars);

%% Create model hierarchy

% The lowest level will be a DMD based model.
load('qgdmd.mat');

modelode{1}  = csl.odetestproblems.qgso.presets.GC('large', 'high');
modelode{1}.TimeSpan = [0 5];
solvermodel{1} = @(~, t, y) deal(t, [y, dmdprop(y)].');
%hm = 0.125;
%solvermodel{1} = @(f, t, y) csl.utils.rk4(f, t, y, round(diff(t)/hm));
%solvermodel{1} = @(f, t, y) ode45(f, t, y);
model{1}  = csl.datools.Model(modelode{1},  solvermodel{1});


% The highest level will be a standard QG model

modelode{2}  = csl.odetestproblems.qgso.presets.GC('large', 'high');
modelode{2}.TimeSpan = [0 5];
solvermodel{2} = @(~, t, y) deal(t, [y, dmdprop(y)].');
hm = 0.125;
solvermodel{2} = @(f, t, y) csl.utils.rk4(f, t, y, round(diff(t)/hm));
%solvermodel{2} = @(f, t, y) ode45(f, t, y);
model{2}  = csl.datools.Model(modelode{2},  solvermodel{2});


%% Create observations

no = 300;

observeindicies = round(linspace(1, model{2}.NumVars, no));
nobsvars = numel(observeindicies);

obserrormodel = csl.datools.error.Gaussian('CovarianceSqrt', speye(nobsvars));
observation = csl.datools.observation.Indexed(model{2}.NumVars, ...
    'ErrorModel', obserrormodel, ...
    'Indicies', observeindicies);
R = speye(nobsvars);

%% Do the rest

ensN = 25;
distfn = model{end}.DistanceFunction;

% We make the assumption that there is no model error
modelerror = csl.datools.error.Error;

load('qgsolargeensemble.mat');
ensembleGenerator = @(x) ysamples(:, randperm(size(ysamples, 2), x));