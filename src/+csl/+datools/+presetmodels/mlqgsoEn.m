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

% create 127 to 63 operator
nf = (nf - 1)/2;
nc = (nf - 1)/2;
rc = 1:nc;
if2cs = [rc, rc, rc];
jf2cs = [(2*rc)-1, 2*rc, (2*rc)+1];
vf2cs = [1/4*ones(1, nc), 1/2*ones(1, nc), 1/4*ones(1, nc)];
If2c1D = sparse(if2cs, jf2cs, vf2cs, nc, nf);
If127_2_c63 = kron(speye(nc), If2c1D) * kron(If2c1D, speye(nf));
Ic63_2_f127 = 4*If127_2_c63.';

nf = (nf - 1)/2;
nc = (nf - 1)/2;
rc = 1:nc;
if2cs = [rc, rc, rc];
jf2cs = [(2*rc)-1, 2*rc, (2*rc)+1];
vf2cs = [1/4*ones(1, nc), 1/2*ones(1, nc), 1/4*ones(1, nc)];
If2c1D = sparse(if2cs, jf2cs, vf2cs, nc, nf);
If63_2_c31 = kron(speye(nc), If2c1D) * kron(If2c1D, speye(nf));
Ic31_2_f63 = 4*If63_2_c31.';

nf = (nf - 1)/2;
nc = (nf - 1)/2;
rc = 1:nc;
if2cs = [rc, rc, rc];
jf2cs = [(2*rc)-1, 2*rc, (2*rc)+1];
vf2cs = [1/4*ones(1, nc), 1/2*ones(1, nc), 1/4*ones(1, nc)];
If2c1D = sparse(if2cs, jf2cs, vf2cs, nc, nf);
If31_2_c15 = kron(speye(nc), If2c1D) * kron(If2c1D, speye(nf));
Ic15_2_f31 = 4*If31_2_c15.';



naturetomodel = csl.datools.observation.Linear(nature.NumVars, 'H', If2c);

%naturetomodel = csl.datools.observation.Observation(nature.NumVars);

%% Create model hierarchy

% The lowest level will be a DMD based model.
load('qgdmd.mat');
load('qgresnn_best.mat');

modelode{1}  = csl.odetestproblems.qgso.presets.GC('large', 'high');
%modelode{1}.Parameters.linearsolver = 'multigrid';
%modelode{1}.Parameters.linearsolvertol = 1e-2;

modelode{1}.TimeSpan = [0 5];
%hm = 0.25;
%solvermodel{1} = @(f, t, y) csl.utils.rk4(f, t, y, round(diff(t)/hm));
%solvermodel{1} = @(f, t, y) ode45(f, t, y);
%solvermodel{1} = @(~, t, y) ode45(@(t, yy) dmdprop{4}(yy), t, y);
%solvermodel{1} = @(~, t, y) deal(t, [y, Ic63_2_f127*dmdprop{5}(If127_2_c63*y)].');
solvermodel{1} = @(~, t, y) deal(t, [y, Ic63_2_f127*(Ic31_2_f63*(Ic15_2_f31*qgresmodel.run(If31_2_c15*(If63_2_c31*(If127_2_c63*y)))))].');
%hm = 0.25;
%solvermodel{1} = @(f, t, y) csl.utils.rk4(@(t, yy) dmdprop{1}(yy), t, y, round(diff(t)/hm));
model{1}  = csl.datools.Model(modelode{1},  solvermodel{1});


% The highest level will be a standard QG model

modelode{2}  = csl.odetestproblems.qgso.presets.GC('large', 'high');
modelode{2}.Parameters.linearsolver = 'multigrid';
modelode{2}.Parameters.linearsolvertol = 1e-4;

modelode{2}.TimeSpan = [0 5];
hm = 0.5;
solvermodel{2} = @(f, t, y) csl.utils.rk4(f, t, y, round(diff(t)/hm));
%solvermodel{2} = @(f, t, y) ode45(f, t, y);
%solvermodel{2} = @(~, t, y) deal(t, [y, dmdprop{3}(y)].');
%solvermodel{2} = @(~, t, y) ode45(@(t, yy) dmdprop{4}(yy), t, y);
%solvermodel{2} = @(~, t, y) deal(t, [y, Ic63_2_f127*dmdprop{5}(If127_2_c63*y)].');
model{2}  = csl.datools.Model(modelode{2},  solvermodel{2});

% Control model for DEnKF
modelodeC  = csl.odetestproblems.qgso.presets.GC('large', 'high');
modelodeC.Parameters = modelode{2}.Parameters;
%modelodeC.Parameters = modelode{1}.Parameters;


modelodeC.TimeSpan = [0 5];
solvermodelC = solvermodel{2};
%solvermodelC = solvermodel{1};
modelC  = csl.datools.Model(modelodeC,  solvermodelC);

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