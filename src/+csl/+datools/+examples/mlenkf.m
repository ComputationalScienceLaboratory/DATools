clear; close all;
% Set rng for standard experiments
rng(17);

addpath('../../ode-test-problems/src');

csl.datools.presetmodels.mlqgsoEn

%inflation = 1.08;
inflation = sqrt(1 + 2/(ensN - 3));

%inflation = 1;

radius = 10;
localization{2} = @(t, y, H) csl.datools.tapering.gc(t, y, radius, distfn, H);
%radius = 20;
%localization{2} = @(t, y, H) csl.datools.tapering.gauss(t, y, radius, distfn, H);

radiuslow = 10;
localization{1} = @(t, y, H) csl.datools.tapering.gc(t, y, radiuslow, distfn, H);
%localization{1} = @(t, y, H) csl.datools.tapering.gauss(t, y, radiuslow, distfn, H);




%model = {model{end}};

ssmall = 1/3;

enkf = csl.datools.statistical.ensemble.MLEnKF(model, ...
    'Observation', observation, ...
    'NumEnsemble', ensN, ...
    'ModelError', modelerror, ...
    'EnsembleGenerator', ensembleGenerator, ...
    'Inflation', inflation, ...
    'Localization', localization, ...
    'Parallel', true, ...
    'Ssmall', ssmall, ...
    'RIPIterations', 0);



enkfC = csl.datools.statistical.ensemble.POEnKF(modelC, ...
    'Observation', observation, ...
    'NumEnsemble', ensN, ...
    'ModelError', modelerror{end}, ...
    'EnsembleGenerator', ensembleGenerator, ...
    'Inflation', inflation, ...
    'Localization', localization{end}, ...
    'Parallel', true, ...
    'RIPIterations', 0);

enkfC.Ensemble = enkf.Ensembles{end};

spinup = 250;
times = 11*spinup;

%spinup = 50;
%times = 300;



msesML = zeros(times - spinup, 1);
msesC = zeros(times - spinup, 1);

rmseML = nan;
rmseC = nan;

ps = '';

% naturefigure = figure;
controlfigure = figure;
% mlfigure = figure;
% topfigure = figure;
botfigure = figure;

for i = 1:times
    
    % forecast
    nature.evolve();
    enkf.forecast();
    enkfC.forecast();

    % observe
    xt = naturetomodel.observeWithoutError(nature.TimeSpan(1), nature.State);
    y = enkf.Observation.observeWithError(model{1}.TimeSpan(1), xt);
    
    % analysis
    enkf.analysis(R, y);
    enkfC.analysis(R, y);
    
    xaTL = mean(enkf.Ensembles{end}, 2);
    
    xaBL = mean(enkf.Ensembles{1}, 2);
    
    xa = enkf.BestEstimate;
    xaC = enkfC.BestEstimate;
    
    
%     set(0,'CurrentFigure', naturefigure)
%     natureode.OutputFcn(nature.TimeSpan(end), xt, []);
%     suptitle('Nature');
%     drawnow;
%     
    set(0,'CurrentFigure', controlfigure)
    natureode.OutputFcn(nature.TimeSpan(end), xaC, []);
    suptitle('Control');
    drawnow;
%     
%     
%     set(0,'CurrentFigure', mlfigure)
%     natureode.OutputFcn(nature.TimeSpan(end), xa, []);
%     suptitle('MLEnKF');
%     drawnow;
%     
%     set(0,'CurrentFigure', topfigure)
%     natureode.OutputFcn(nature.TimeSpan(end), xaTL, []);
%     suptitle('Top Layer');
%     drawnow;
    
    set(0,'CurrentFigure', botfigure)
    natureode.OutputFcn(nature.TimeSpan(end), xaBL, []);
    suptitle('Bottom Level');
    drawnow;
    
    
    %subplot(2, 2, 1);
    %imagesc(reshape(nature.State, 255, 255));
    %imagesc(reshape(xt, 127, 127).');
    %axis square; colorbar;
    %title('Nature');
    %     subplot(2, 2, 2);
    %     imagesc(reshape(xaC, 127, 127).');
    %     axis square; colorbar;
    %     title('EnKF');
    %
    %
    %     subplot(2, 3, 4);
    %     imagesc(reshape(xa, 127, 127).');
    %     axis square; colorbar;
    %     title('MLEnKF');
    %
    %     subplot(2, 3, 5);
    %     imagesc(reshape(xaTL, 127, 127).');
    %     axis square; colorbar;
    %     title('MLEnKF (Top Layer)');
    %
    %     subplot(2, 3, 6);
    %     imagesc(reshape(xaBL, 127, 127).');
    %     axis square; colorbar;
    %     title('MLEnKF (Bottom Layer)');
    %     drawnow;

    
    if i > spinup
        msesML(i - spinup) = mean((xaTL - xt).^2);
        rmseML = sqrt(mean(msesML(1:(i - spinup))));
        
        msesC(i - spinup) = mean((xaC - xt).^2);
        rmseC = sqrt(mean(msesC(1:(i - spinup))));
    end
    
    for kk = 1:numel(ps)
        fprintf('\b');
    end
    
    ps = sprintf('step: %d, rmseML: %.5f, rmseC: %.5f\n', i, rmseML, rmseC);
    
    fprintf(ps);
end
