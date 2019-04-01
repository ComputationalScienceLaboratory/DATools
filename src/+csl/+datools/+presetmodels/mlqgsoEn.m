
DT = 5;


%% Create Nature

load('qglargeensemble.mat');

natureode = otp.qg.presets.PopovSandu2019('Nature');
%natureode = csl.odetestproblems.qgso.presets.GC('large', 'high');

natureode.Y0 = ysamples(:, end);

natureode.TimeSpan = [0, DT];



solvernature = @(f, t, y) ode45(f, t, y);

hm = 0.25;
%solvernature = @(f, t, y) csl.utils.rk4(f, t, y, round(diff(t)/hm));
solvernature = @(f, t, y) csl.utils.tvdrk3(f, t, y, round(diff(t)/hm));


% nature will have no sythetic error
naturesyntherror = csl.datools.error.Error;
nature = csl.datools.Model(natureode, solvernature, ...
    'SynthError', naturesyntherror);

naturetomodel = csl.datools.observation.Observation(nature.NumVars);

%% Create model hierarchy

% The lowest level will be a DMD based model.
%load('qgresnn_best.mat');
load('qgdmd.mat');
load('qg_pod_basis.mat');
load('qgsqrtcov.mat')


k = 90;
U = U(:, 1:k);
UArakawa = UArakawa(:, 1:k);
PArakawa = PArakawa(:, 1:k);

UDMD = UDMD(:, 1:k);
SDMD = SDMD(1:k, 1:k);

modelode{1} = otp.qg.presets.PopovSandu2019('POD');



modelode{1}.TimeSpan = [0 DT];

% POD
modelode{1}.Parameters.pod = struct('basis', U, 'arakawabasis', UArakawa, 'arakawaDEIM', PArakawa);
hm = 0.25;
solvermodel{1} = @(~, t, y) modelode{1}.fromPOD(csl.utils.rk4_filter(modelode{1}.RhsPOD.F, t, U'*y, round(diff(t)/hm), modelode{1}.FilterPOD));
solvermodel{1} = @(~, t, y) modelode{1}.fromPOD(csl.utils.tvdrk3_filter(modelode{1}.RhsPOD.F, t, U'*y, round(diff(t)/hm), modelode{1}.FilterPOD));


%solvermodel{1} = @(~, t, y) modelode{1}.fromPOD(csl.utils.rk4_filter(@(tt, yy) modelode{1}.RhsPOD.F(tt, yy) + U'*PODbias/DT, t, U'*y, round(diff(t)/hm), modelode{1}.FilterPOD));


%solvermodel{1} = @(~, t, y) modelode{1}.fromPOD(csl.utils.rk4(modelode{1}.RhsPOD.F, t, U'*y, round(diff(t)/hm)));



%DMD
%modelode{1}.Parameters.pod = struct('basis', UDMD, 'arakawabasis', UArakawa, 'arakawaDEIM', PArakawa);
%hm = 0.25;
%solvermodel{1} = @(f, t, y) modelode{1}.fromPOD(csl.utils.rk4_filter(@(t, yy) SDMD*yy, t, UDMD'*y, round(diff(t))/hm, modelode{1}.FilterPOD));

%solvermodel{1} = @(f, t, y) csl.utils.rk4_filter(@(t, yy) dmdprop{end}(yy), t, y, round(diff(t))/hm, modelode{1}.Filter);


model{1}  = csl.datools.Model(modelode{1},  solvermodel{1});


% The highest level will be a standard QG model

modelode{2} = otp.qg.presets.PopovSandu2019('Model');

%modelode{2}.Parameters.e = 1e-7;
%modelode{2}.Parameters.linearsolver = 'multigrid';
%modelode{2}.Parameters.linearsolvertol = 1e-7;

modelode{2}.TimeSpan = [0 DT];
hm = 0.25;
solvermodel{2} = @(f, t, y) csl.utils.rk4(f, t, y, round(diff(t)/hm));
solvermodel{2} = @(f, t, y) csl.utils.tvdrk3(f, t, y, round(diff(t)/hm));
%solvermodel{2} = @(f, t, y) ode45(f, t, y);
%solvermodel{2} = @(~, t, y) deal(t, [y, dmdprop{3}(y)].');
%solvermodel{2} = @(~, t, y) ode45(@(t, yy) dmdprop{4}(yy), t, y);
%solvermodel{2} = @(~, t, y) deal(t, [y, Ic63_2_f127*dmdprop{5}(If127_2_c63*y)].');
model{2}  = csl.datools.Model(modelode{2},  solvermodel{2});

% Control model for DEnKF
modelodeC = otp.qg.presets.PopovSandu2019('Model');
modelodeC.Parameters = modelode{2}.Parameters;
%modelodeC.Parameters = modelode{1}.Parameters;

modelodeC.TimeSpan = [0 DT];
solvermodelC = solvermodel{2};
%solvermodelC = solvermodel{1};
modelC  = csl.datools.Model(modelodeC,  solvermodelC);


% MAGIC REMOVE
%model = fliplr(model);


%% Create observations

no = 150;
%no = 50;

% we only want to observe the top layer
nvar = model{2}.NumVars/model{2}.ODEModel.Parameters.nlayers;

observeindicies = round(linspace(1, nvar, no));
nobsvars = numel(observeindicies);

R = 4*speye(nobsvars);

obserrormodel = csl.datools.error.Gaussian('CovarianceSqrt', sqrtm(R));
observation = csl.datools.observation.Indexed(model{2}.NumVars, ...
    'ErrorModel', obserrormodel, ...
    'Indicies', observeindicies);


indexing = csl.datools.observation.Indexed(model{2}.NumVars, ...
    'Indicies', observeindicies);
H = indexing.linearization([], []);

flowvel = @(~, psi) H*natureode.FlowVelocityMagnitude(psi);
linflowvel = @(~, psi) (H*natureode.JacobianFlowVelocityMagnitude(psi));

% observation = csl.datools.observation.Nonlinear(model{2}.NumVars, ...
%      'ErrorModel', obserrormodel, ...
%      'F', flowvel, ...
%      'J', linflowvel);


%% Do the rest

ensN = 25;
distfn = model{end}.DistanceFunction;

% We make the assumption that there is no model error
modelerror{2} = csl.datools.error.Error;

% However for the POD model there is model error
%modelerror{1} = csl.datools.error.Gaussian('CovarianceSqrt', Qsqrt);
modelerror{1} = csl.datools.error.Error;
%modelerror{1} = csl.datools.error.Gaussian('CovarianceSqrt', Qsqrt, 'Bias', PODbias);
%modelerror{1} = csl.datools.error.Gaussian('CovarianceSqrt', 0*Qsqrt, 'Bias', PODbias);




ensembleGenerator = @(x) ysamples(:, randperm(size(ysamples, 2), x));