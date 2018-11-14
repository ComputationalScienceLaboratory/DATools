function LF = pa(rmin, rmax, nrs, rsfun, distfun, m, locfun)

LF = @(obj) pamat(obj, locfun, rsfun, distfun, obj.Problem.TimeSpan(1), obj.Problem.Y0, m, rmin, rmax, nrs);

end

function rho = pamat(obj, locfun, rsfun, distfun, t, y, m, rmin, rmax, nrs)

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

Jn = @(r) costfun(Pf, locfun, numVars, rsfun(r), distfun, t, y, m, d, H, R, Xf, ensN);


% rss = rmin:rmax;
% 
% Jss = arrayfun(Jn, rss);
% 
% %plot(rss, Jss);
% %drawnow;
% 
% [~, ri] = min(Jss);
% 
% r = rss(ri);

psoptions = optimoptions('patternsearch', 'Display', 'off');
r = patternsearch(Jn, (rmins + rmaxs)/2, [], [], [], [], rmins, rmaxs, [], psoptions);



rs = rsfun(r);
rho = locfun(numVars, rs, distfun, t, y, m);

end

function J = costfun(Pf, locfun, n, rs, distfun, t, y, m, d, H, R, Xf, ensN)

Pfr = (locfun(n, rs, distfun, t, y, m) .* Pf);


F = Pfr(H, H);
S = (F + R);

K1 = Pf(:, H)/S;

Xa = Xf + K1 * d;

Xam = mean(Xa, 2);

tmp = Xa - Xam * ones(1, ensN);
Pa1 = (locfun(n, rs, distfun, t, y, m) .* (1/(ensN - 1) *(tmp*tmp')));

K2 = Pfr(:, H)/S;

Pa2 = Pfr - K2 * Pfr(H, :);


% subplot(1,2,1)
% imagesc(Pa1);
% subplot(1,2,2)
% imagesc(Pa2);
% 
% drawnow

dd = Pa1 - Pa2;

M = length(dd);

J = norm(dd);

end

