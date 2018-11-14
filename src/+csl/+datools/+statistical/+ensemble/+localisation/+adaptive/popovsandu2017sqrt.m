function LF = popovsandu2017sqrt(rmin, rmax, nrs, rsfun, distfun, m, locfun)

LF = @(obj) popovsandu2017mat(obj, B, locfun, rsfun, distfun, obj.Problem.TimeSpan(1), obj.Problem.Y0, m, rmin, rmax, nrs);

end

function rs = popovsandu2017mat(obj, B, locfun, rsfun, distfun, t, y, m, rmin, rmax, nrs)

ensN = obj.EnsembleSize;
Xf =  obj.CurrentBestGuess;
Xfm = mean(Xf, 2);
R = obj.CurrentObsErr;
H = obj.H;

numVars = obj.Problem.NumVars;
y = obj.CurrentObs;

if mod(numVars, nrs) ~= 0
    error('Invalid value for number of radii');
end

rmins = rmin * ones(nrs, 1);
rmaxs = rmax * ones(nrs, 1);

Jn = @(r) costfun(obj, B, locfun, numVars, rsfun(r), distfun, t, Xf, Xfm, y, H, R);

%psoptions = optimoptions('patternsearch', 'Display', 'off', 'MaxTime', 0.25, ...
%    'UseCompletePoll', true, ...
%    'UseParallel', true);
psoptions = optimoptions('patternsearch', 'Display', 'off', 'MaxTime', 2.5);
r = patternsearch(Jn, (rmins + rmaxs)/2, [], [], [], [], rmins, rmaxs, [], psoptions);

% fminconoptions = optimoptions('fmincon', 'Display', 'off');
% r = fmincon(Jn, (rmins + rmaxs)/2, [], [], [], [], rmins, rmaxs, [], fminconoptions);


rs = rsfun(r);

end

function J = costfun(obj, B, locfun, n, rs, distfun, t, Xf, Xfm, y, H, R)

% Reich and Cotter
Af = Xf - Xfm*ones(1, ensN);

HAf = Af(H, :);

Xa = zeros(size(Xf));

for k = 1:numel(Xfm)
    %
    
    Ctilde = generateCTilde(obj, k);
    
    HCtildeH = Ctilde(H, H);
    
    HCtildeHRi = (HCtildeH/R);
    
    Sk2i = eye(ensN, ensN) + 1/(ensN - 1) * (HAf.') * (HCtildeHRi) * HAf;
    
    Ski = sqrtm(Sk2i);
    
    Aak = Af(k, :)/Ski;
    
    wk = 1/ensN - 1/(ensN - 1) * (Sk2i\(HAf.')) * HCtildeHRi * (Xfm(H, :) - y);
    
    
    Xakm = Xf(k, :)*wk;
    
    
    Xa(k, :) = Xakm +Aak;
end



XaminusXf = Xa - Xf;
yminusHXa = y - Xa(H, :);

Jf = 0.5*norm((XaminusXf.')/B*(XaminusXf));
Jo = 0.5*norm((yminusHXa.')/R*(yminusHXa));

J = Jf + Jo;

end

function Ctilde = generateCTilde(obj, k)

Ctilde = obj.LocalisationFunction(obj, k);

end
