clear;


%% Model definition
model = csl.odetestproblems.qgso.presets.GC('large', 'high');

%% First model evolution
[t, y] = ode45(model.F, [0 1000], model.Y0);
model.Y0 = y(end, :).';

%% Second model evolution
tstart = t(end);
h = 0.1;
tend = t(end) + 1000;
[~, y] = ode45(model.F, tstart:h:tend, model.Y0);

%% Naive linear DMD algorithm
VNm1 = y(1:(end - 1), :).';
VN   = y(2:end, :).';
dVN = (VN - VNm1)/h;
ks = [100 200 500 1000];
dmdprop = cell(numel(ks), 1);
for ki = 1:numel(ks)
    
    k = ks(ki);
    [U, Sigma, W] = svds(VNm1, k);
    
    %Stilde = U'*VN*W/Sigma;
    Stilde = U'*dVN*W/Sigma;
    
    dmdprop{ki} = @(v) U*(Stilde*(U'*v));
    
    ds = dmdprop{ki}(VNm1) - dVN;
    
    rms(ds(:))
    
end

%% Variable saving

save('qgdmd.mat', 'dmdprop', 'ks');
