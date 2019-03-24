%% Create Nature

natureode = otp.qg.presets.PopovSandu2019('Nature');
%natureode = csl.odetestproblems.qgso.presets.GC('large', 'high');

natureode.TimeSpan = [0, 5];

solvernature = @(f, t, y) ode45(f, t, y);

% nature will have no sythetic error
naturesyntherror = csl.datools.error.Error;
nature = csl.datools.Model(natureode, solvernature, ...
    'SynthError', naturesyntherror);

naturetomodel = csl.datools.observation.Observation(nature.NumVars);

%% Create model hierarchy

% The lowest level will be a DMD based model.
%load('qgresnn_best.mat');
load('qg_pod_basis.mat');

modelode{1} = otp.qg.presets.PopovSandu2019('POD');



modelode{1}.TimeSpan = [0 5];

% POD
modelode{1}.Parameters.pod = struct('basis', U, 'arakawabasis', UArakawa, 'arakawaDEIM', PArakawa);

%solvermodel{1} = @(~, t, y) modelode{1}.fromPOD(ode45(modelode{1}.RhsPOD.F, t, U'*y));
hm = 1.25;
solvermodel{1} = @(~, t, y) modelode{1}.fromPOD(csl.utils.rk4_filter(modelode{1}.RhsPOD.F, t, U'*y, round(diff(t)/hm), modelode{1}.Filter));

%solvermodel{1} = @(~, t, y) modelode{1}.fromPOD(csl.utils.rk4(modelode{1}.RhsPOD.F, t, U'*y, round(diff(t)/hm)));


%hm = 0.25;
%solvermodel{1} = @(f, t, y) csl.utils.rk4(@(t, yy) dmdprop{1}(yy), t, y, round(diff(t)/hm));
model{1}  = csl.datools.Model(modelode{1},  solvermodel{1});


% The highest level will be a standard QG model

modelode{2} = otp.qg.presets.PopovSandu2019('Model');

%modelode{2}.Parameters.e = 1e-7;
%modelode{2}.Parameters.linearsolver = 'multigrid';
%modelode{2}.Parameters.linearsolvertol = 1e-7;

modelode{2}.TimeSpan = [0 5];
hm = 0.25;
solvermodel{2} = @(f, t, y) csl.utils.rk4(f, t, y, round(diff(t)/hm));
%solvermodel{2} = @(f, t, y) ode45(f, t, y);
%solvermodel{2} = @(~, t, y) deal(t, [y, dmdprop{3}(y)].');
%solvermodel{2} = @(~, t, y) ode45(@(t, yy) dmdprop{4}(yy), t, y);
%solvermodel{2} = @(~, t, y) deal(t, [y, Ic63_2_f127*dmdprop{5}(If127_2_c63*y)].');
model{2}  = csl.datools.Model(modelode{2},  solvermodel{2});

% Control model for DEnKF
modelodeC = otp.qg.presets.PopovSandu2019('Model');
modelodeC.Parameters = modelode{2}.Parameters;
%modelodeC.Parameters = modelode{1}.Parameters;

modelodeC.TimeSpan = [0 5];
solvermodelC = solvermodel{2};
%solvermodelC = solvermodel{1};
modelC  = csl.datools.Model(modelodeC,  solvermodelC);

%% Create observations

% no = 300;
no = 300;

% we only want to observe the top layer
nvar = model{2}.NumVars/model{2}.ODEModel.Parameters.nlayers;

observeindicies = round(linspace(1, model{2}.NumVars, no));
nobsvars = numel(observeindicies);

R = 1*speye(nobsvars);

obserrormodel = csl.datools.error.Gaussian('CovarianceSqrt', sqrtm(R));
observation = csl.datools.observation.Indexed(model{2}.NumVars, ...
    'ErrorModel', obserrormodel, ...
    'Indicies', observeindicies);


%% Do the rest

ensN = 12;
distfn = model{end}.DistanceFunction;

% We make the assumption that there is no model error
modelerror = csl.datools.error.Error;

load('qglargeensemble.mat');
ensembleGenerator = @(x) ysamples(:, randperm(size(ysamples, 2), x));