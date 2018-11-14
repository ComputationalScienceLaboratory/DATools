function LF = popovsandu2017(rmin, rmax, nrs, rsfun, distfun, m, locfun)

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

%psoptions = optimoptions('patternsearch', 'Display', 'off', 'MaxTime', 0.25, ...
%    'UseCompletePoll', true, ...
%    'UseParallel', true);
%psoptions = optimoptions('patternsearch', 'Display', 'off', 'MaxTime', 2.5);
%r = patternsearch(Jn, (rmins + rmaxs)/2, [], [], [], [], rmins, rmaxs, [], psoptions);


%% ALTERNATIVE

r = (rmins + rmaxs)/2;

passes = 2;

% plot(1:numel(r), r, '-o');
% axis([0 numel(r)+1 rmin-1 rmax+1]);
% title('Pass=1, ri=1');
% drawnow;

rsperm = randperm(numel(r));

for pass = 1:passes
    for ii = rsperm
        replace = @(xs, jj, x) [xs(1:(jj-1)); x; xs((jj+1):end)];
        Jnn = @(nr) Jn(replace(r, ii, nr));
        allrs = rmin:0.25:rmax;
        [~, I] = min(arrayfun(Jnn, allrs));
        r = replace(r, ii, allrs(I)); 
%         plot(1:numel(r), r, '-o');
%         axis([0 numel(r)+1 rmin-1 rmax+1]);
%         title(sprintf('Pass=%d, ri=%d', pass, ii));
%         drawnow;
    end
end


%% END


%% Grad Desc

% alpha = 1;
% 
% r = (rmins + rmaxs)/2;
% 
% k = 1;
% pp = 2;
% q = 0;
% 
% epsilon = sqrt(eps);
% 
% for its = 1:10
% 
% Jnr = Jn(r);
% % A = @(b) (Jn(r + epsilon*b) - Jnr)/epsilon;
% 
% % [U, S, V] = csl.dataassimilation.utils.martSVDmatfree(A, nrs, k, pp, q);
% 
% % g = U*S*V';
% 
% g = zeros(nrs, 1);
% 
% for i = 1:nrs
%     e = zeros(nrs, 1);
%     e(i) = 1;
%     g(i) = (Jn(r + epsilon*e) - Jnr)/epsilon;
% end
% p = -g;
% 
% Phi0 = Jnr;
% Phi1 = Jn(r + p);
% dPhi0 = p.' * g;
% 
% alpha = -dPhi0/(2*(Phi1 - Phi0 - dPhi0));
% 
% alpha = max(alpha, 0.1);
% 
% r = r + alpha*p;
% 
% end
% 
% r = max(r, 0.5);
% r = min(r, 20);
% 
% plot(1:numel(r), r, '-o');
% axis([0 numel(r)+1 0.5 20]);
% drawnow;




%% END


% fminconoptions = optimoptions('fmincon', 'Display', 'off');
% r = fmincon(Jn, (rmins + rmaxs)/2, [], [], [], [], rmins, rmaxs, [], fminconoptions);


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

