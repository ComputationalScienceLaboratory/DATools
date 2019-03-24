clear;
% Set rng for standard experiments
rng(20);

addpath('../../ode-test-problems/src');

csl.datools.presetmodels.mlqgsoEn


xt = naturetomodel.observeWithoutError(nature.TimeSpan(1), nature.State);

clf;

xt = (If63_2_c31*(If127_2_c63*xt));




for i = 1:1000
    %imagesc(reshape(xt, 127, 127));
    imagesc(reshape(xt, 31, 31));
    axis square; colorbar;
    title('Model');
    drawnow;
    
    
    %xt = If31_2_c15*(If63_2_c31*(If127_2_c63*xt));
    
    %xt = qgresmodel.run(xt);

    %xt = Ic63_2_f127*(Ic31_2_f63*(Ic15_2_f31*xt));
    
    
    
    
    %xt = (If63_2_c31*(If127_2_c63*xt));
    
    xt = qgresmodel.run(xt);

    %xt = Ic63_2_f127*(Ic31_2_f63*(xt));
    
end
