function LF = popovsandu20174D(rmin, rmax, nrs, rsfun, distfun, m, locfun)

LF = @(obj) popovsandu2017mat(obj, locfun, rsfun, distfun, obj.Problem.TimeSpan(1), obj.Problem.Y0, m, rmin, rmax, nrs);

end

function rho = popovsandu2017mat(obj, locfun, rsfun, distfun, t, y, m, rmin, rmax, nrs)

ensN = obj.EnsembleSize;
Xf =  obj.CurrentBestGuess;
R = obj.CurrentObsErr;
H = obj.H;

numVars = obj.Problem.NumVars;
Ys = obj.CurrentSyntheticObs;

if mod(numVars, nrs) ~= 0
    error('Invalid value for number of radii');
end

Xfm = mean(Xf, 2);

tmp = Xf - Xfm*ones(1,ensN);
Pf = 1/(ensN - 1) *(tmp*tmp');

rmins = rmin * ones(nrs, 1);
rmaxs = rmax * ones(nrs, 1);

d = Ys - Xf(H, :);

Jn = @(r) costfun(Pf, locfun, numVars, rsfun(r), distfun, t, y, m, d, H, R);

psoptions = optimoptions('patternsearch', 'Display', 'off');
r = patternsearch(Jn, (rmins + rmaxs)/2, [], [], [], [], rmins, rmaxs, [], psoptions);

rs = rsfun(r);
rho = locfun(numVars, rs, distfun, t, y, m);

% Pfr = (rho .* Pf);
% disp('------');
% M = length(Pfr);
% (det(Pfr))^(-(M+1)/2)
% (det(R))^(-(M+1)/2)

end

function J = costfun(Pf, locfun, n, rs, distfun, t, y, m, d, H, R)

Pfr = (locfun(n, rs, distfun, t, y, m) .* Pf);

F = Pfr(H, H);
S = (F + R);

W = F/S;
G = d - W*d;

Jf = 1/2*norm(d'*(S\W)*d);
Jo = 1/2*norm(G'/R*G);

% Jf = sqrt(det(Pfr)) * 1/2*norm(d'*(S\W)*d);
% Jo = sqrt(det(R))   * 1/2*norm(G'/R*G);

% M = length(Pfr);
% Jf = (det(Pfr))^((M+1)/2) * 1/2*norm(d'*(S\W)*d);
% Jo = (det(R))^((M+1)/2)   * 1/2*norm(G'/R*G);

J = Jf + Jo;

end

