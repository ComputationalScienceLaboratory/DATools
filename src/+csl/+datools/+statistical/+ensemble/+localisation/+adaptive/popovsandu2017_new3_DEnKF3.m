function LF = popovsandu2017_new3_DEnKF2(rmin, rmax, nrs, rsfun, distfun, m, locfun, betamean, D, inflation)

LF = @(obj, H) popovsandu2017mat(obj, locfun, rsfun, distfun, obj.Problem.TimeSpan(1), obj.Problem.Y0, m, rmin, rmax, nrs, betamean, D, H, inflation);

end

function rhoHt = popovsandu2017mat(obj, locfun, rsfun, distfun, t, y, m, rmin, rmax, nrs, betamean, D, H, inflation)

ensN = obj.EnsembleSize;
Xf =  obj.CurrentBestGuess;
R = obj.CurrentObsErr;
% H = obj.H;

numVars = obj.Problem.NumVars;
%Ys = obj.CurrentSyntheticObs;
Y = obj.CurrentObs;

if mod(numVars, nrs) ~= 0
    error('Invalid value for number of radii');
end

Xfm = mean(Xf, 2);

Af = Xf - Xfm*ones(1, ensN);
% Af = inflation * Af;
% Xf = Xfm*ones(1, ensN) + Af;

HPfHt = 1/(ensN - 1) *(Af(H, :)*(Af(H, :)'));

%rmins = rmin * ones(nrs, 1);
%rmaxs = rmax * ones(nrs, 1);

d = Y - Xfm(H, :);



% betafun = @(b) (rmax - rmin)*(erf(b) + 1)/2 + rmin;
% betafuninv = @(bi) erfinv(2*(bi - rmin)/(rmax - rmin) - 1);

betafun = @(b) b;
betafuninv = @(bi) bi;

%Jn = @(b) costfun(Pf, locfun, numVars, b, rsfun, distfun, t, y, m, d, H, R, betamean, D, Xf, Y);

% psoptions = optimoptions('patternsearch', 'Display', 'off', 'MaxTime', 2.5, ...
%    'UseCompletePoll', true, ...
%    'UseParallel', true);

% options = optimoptions('fmincon', 'Display', 'off', ...
%    'UseParallel', true);
% beta = fmincon(Jn, betamean * ones(nrs, 1), [], [], [], [], rmins, rmaxs, [], options);


% psoptions = optimoptions('patternsearch', 'Display', 'off', 'MaxTime', 1);
% beta = patternsearch(Jn, betamean * ones(nrs, 1), [], [], [], [], rmins, rmaxs, [], psoptions);


% psoptions = optimoptions('fminunc', 'Display', 'off', ...
%    'UseParallel', true);
% beta = fminunc(Jn, betamean * ones(nrs, 1), psoptions);


% psoptions = optimoptions('particleswarm', ...
%     'Display', 'off', ...
%     'UseParallel', true, ...
%     'MaxIterations', 10);
% beta = particleswarm(Jn, nrs, rmins, rmaxs, psoptions);




%% ALTERNATIVE

%
% beta = betamean * ones(nrs, 1);
% if nrs == 1
%     passes = 1;
% else
%     passes = 1;
% end
%
% inflation = 1;
%
% Jn = @(b) costfun(Pf, locfun, numVars, b, inflation, rsfun, distfun, t, y, m, d, H, R, betamean, D, Xf, Y);
%
% rsperm = randperm(nrs);
%
% for pass = 1:passes
%     for ii = rsperm
%         replace = @(xs, jj, x) [xs(1:(jj-1)); x; xs((jj+1):end)];
%         Jnn = @(nb) Jn(replace(beta, ii, nb));
%         allrs = betafuninv((rmin):0.5:(rmax));
%
%         crs = zeros(numel(allrs, 1));
%
%         parfor allrsi = 1:numel(allrs)
%             crs(allrsi) = Jnn(allrs(allrsi));
%         end
%
%         [~, I] = min(crs);
%
%         %[~, I] = min(arrayfun(Jnn, allrs));
%         beta = replace(beta, ii, allrs(I));
%     end
% end

% with inflation

beta = betamean * ones(nrs, 1);
if nrs == 1
    passes = 1;
else
    passes = 1;
end

Jn = @(b) costfun(HPfHt, locfun, numVars, b, rsfun, distfun, t, y, m, d, H, R, betamean, D, Xf, Y, Af, inflation, rmin, rmax);

% rsperm = randperm(numel(beta));




% psoptions = optimoptions('fminunc', 'Display', 'off', ...
%     'UseParallel', false, 'Algorithm', 'quasi-newton');
% beta = fminunc(Jn, beta, psoptions);

% optoptions = optimoptions('fmincon', 'Display', 'off', ...
%     'UseParallel', false, 'MaxFunctionEvaluations', 40);
% beta = fmincon(Jn, beta, [], [], [] ,[] ,[] ,[], [], optoptions);

% 
optoptions = optimoptions('fmincon', 'Display', 'off', ...
    'UseParallel', false, 'MaxFunctionEvaluations', 5);
rsperm = 1:numel(beta);
for pass = 1:passes
    for ii = rsperm
        %replace = @(xs, jj, x) [xs(1:(jj-1)); x; xs((jj+1):end)];
        Jnn = @(nb) Jn([beta(1:(ii-1)); nb; beta((ii+1):end)]);

        zz = fmincon(Jnn, beta(ii), [], [], [] ,[] ,[] ,[], [], optoptions);
        
        %[~, I] = min(arrayfun(Jnn, allrs));
        beta = [beta(1:(ii-1)); zz; beta((ii+1):end)];
    end
end


% rsperm = 1:numel(beta);
%
% for pass = 1:passes
%     for ii = rsperm
%         replace = @(xs, jj, x) [xs(1:(jj-1)); x; xs((jj+1):end)];
%         Jnn = @(nb) Jn(replace(beta, ii, nb));
%
%         allrs = [betafuninv((rmin):0.5:(rmax))];
%
%         %allrs = [10:5:40];
%
%
%         crs = zeros(numel(allrs, 1));
%
%         parfor allrsi = 1:numel(allrs)
%             crs(allrsi) = Jnn(allrs(allrsi));
%         end
%
%         [~, I] = min(crs);
%         %[~, I] = min(arrayfun(Jnn, allrs));
%         beta = replace(beta, ii, allrs(I));
%     end
% end


%%%%%update
% global betam;
% global Dm;
% global Nm;
%
% if isempty(Nm)
%     betam = 2.5;
%     Dm = 10;
%     Nm = 50;
% end
%
%
% betamean = betam;
% D = Dm;
% delta = beta - betam;
% Dm = 1/Nm * ((Nm - 1)*Dm + delta*delta');
%
% betam = (betam*Nm + beta)/(Nm + 1);
% Nm = Nm + 1;

% global zzz;
% global rads;
% if isempty(zzz)
%     clf;
%     zzz = 1;
%     rads = [beta];
% else
%     zzz = zzz + 1;
%     rads = [rads, beta];
% end
% if mod(zzz, 25) == 0
%     clf;
%     hold all;
%     plot(1:zzz, rads, '.k');
%     plot([0 zzz], [betamean betamean], '-b');
%     axis([0 zzz, 0, 20]);
%     title(sprintf('r = %.5f, bm = %.5f, bmr=%.5f, D = %.5f', beta, betamean, betafun(betamean), D));
%     drawnow;
% end

% bs = linspace(rmin, rmax);
% [Js, Jfs, Jos] = arrayfun(Jn, bs);
% clf;
% hold all;
% plot(bs, Js);
% plot(bs, Jfs);
% plot(bs, Jos);
% pause;

%Jn(beta);


% global zzz2;
% if isempty(zzz2)
%     zzz2 = [];
% end
% if beta ~= rmax
%     zzz2 = [zzz2 beta];
% end
%
% z = tabulate(zzz2);
% clf;
% plot(z(:, 1), z(:, 3)/100);
% if mod(numel(zzz2), 10) == 0
%     drawnow;
% end


%
% global zzz;
% if isempty(zzz)
%     clf;
%     zzz = 1;
% else
%     zzz = zzz + 1;
% end
% hold all;
% plot(zzz, r, '.k');
% plot([0 zzz], [2.5 2.5], '-b');
% axis([0 zzz, 0, 10]);
% title(sprintf('r = %f', r));
% if mod(zzz, 10) == 0
%     drawnow;
% end


%% END


%% Grad Desc

% beta = betamean * ones(nrs, 1);
%
% epsilon = sqrt(eps);
%
% for its = 1:50
%
% Jnr = Jn(beta);
%
% g = zeros(nrs, 1);
%
% for i = 1:nrs
%     e = zeros(nrs, 1);
%     e(i) = 1;
%     g(i) = (Jn(beta + epsilon*e) - Jnr)/epsilon;
% end
% p = -g;
%
% Phi0 = Jnr;
% Phi1 = Jn(beta + p);
% dPhi0 = p.' * g;
%
% alpha = -dPhi0/(2*(Phi1 - Phi0 - dPhi0));
%
% alpha = max(alpha, 0.1);
%
% beta = beta + alpha*p;
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


%rs = rsfun(betafun(beta));

rhoHt = locfun(numVars, beta, distfun, t, y, m, 1, H);
%rhoHt = locfun(numVars, rs, distfun, t, y, m, 1, H);

% title(sprintf('%.5e', beta));
% drawnow;

% Pfr = (rho .* Pf);
% disp('------');
% M = length(Pfr);
% (det(Pfr))^(-(M+1)/2)
% (det(R))^(-(M+1)/2)

end

function [J, Jf, Ju, Jo] = costfun(HPfHt, locfun, n, beta, rsfun, distfun, t, y, m, d, H, R, bm, Di, Xf, Y, Af, inflation, rmin, rmax)
% fprintf('%.5e\n',beta);
% Y = mean(Ys, 2);

F = (locfun(n, rsfun(beta), distfun, t, y, m, 1, H, H) .* HPfHt);

ensN = size(Xf, 2);

%F = PfHtr(H, :);
S = (F + R);

HK = F/S;


HAf = Af(H, :);

gs = (speye(numel(d)) - HK)*(d*ones(1, ensN) - HAf) - (1/2)*HK*HAf;

%gs = d*ones(1, ensN) + (1/2)*HK*(Af(H, :));

%W = K(H, :);

% Xam = Xfm + Kd
% Aa  = Af - 0.5*K*H*Af;



%W = F/S;
%G = d - W*d;

%Jf = 1/2*norm(d'*(S\W)*d);

%Jo = 1/2*norm(G'/R*G);


% Normal
%D = Di*eye(numel(beta));
%Jf = 1/2*norm((beta-bm).'/D*(beta-bm));



%Xfm = mean(Xf, 2);
%K = Pfr(:, H)/S;
%Xam = Xfm + K*d;
%Aa  = Af - 0.5*K*Af(H, :);

%Aa = inflation * Aa;

%Xa = Xam*ones(1, ensN) + Aa;

% Xa - Xf = Xam*1 + Aa - Xfm*1 - Af
%         = K*(y-H*Xfm)*1 - 0.5*K*H*Af
%         = PfHtSi*z
% z = (y-H*Xfm)*1 - 0.5*H*Af

z = d*ones(1, ensN) - 0.5*Af(H, :);


n = size(Xf,1);

% Gamma
%aa = 9.640796014193920e-01;
%bb = 1/5.421435022127940e+00;



%rm = aa/bb;

%rm = bm;

%Pfrm = (locfun(n, rsfun(beta), distfun, t, y, m) .* PfHt);
%
%Krm = Pfrm(:, H)/(Pfrm(H, H) + R);

%Xarm = Xf + Krm*d;
%Xarmm = mean(Xa, 2);



%tmp = Xarm - Xarmm*ones(1, ensN);
%Parm = 1/(ensN - 1) *(tmp*tmp');

%Parm = (locfun(n, rsfun(rm), distfun, t, y, m) .* Parm);

%Xam = mean(Xa, 2);
%Xfm = mean(Xf, 2);

% THIS MIGHT BE GOOD BUT ALT FOR NOW
Jf = 0;
for asd = 1:ensN
    %         de = d(:, asd);
    ze = z(:, asd);
    zee = S\ze;
    ZZ = zee.'*F*zee;
    %ZZ = ze.'*(S\F/S)*ze;
    %XaXfd = Xa(:, asd) - Xf(:, asd);
    %Jft = 1/2*XaXfd.'/Pfrm*XaXfd;
    Jf1 = ZZ;
    
    Jf = Jf + Jf1;
end

Jf = 1/2*(Jf);


%     dm = mean(d, 2);
%
%     ZZ = dm.'*(S\F/S)*dm;
%
%     Jf = 1/2*ZZ;


%Jf = 1/2*norm((Xam - Xfm).'/Parm*(Xam - Xfm));
%Jf = 1/2*norm((Xam - Xfm).'/Pfr*(Xam - Xfm));

gam = rmin - 0.5;
theta = Di/(bm - gam);
k = ((bm - gam)^2)/Di;
aa= k;
bb = 1/theta;

if Di == inf
    if any(beta > rmin)
        Ju = 0;
    else
        Ju = inf;
    end
else
    if any(beta <= gam) || any(beta > rmax)
        Ju = inf;
    else
        Ju = (  (-sum((aa-1)*log(beta - gam)) + sum((beta - gam)*bb))  );
    end
end





%mu = bm;
%lambda = (mu^3)/Di;
%Ju = ( lambda*((beta-mu)^2) )/( 2*(mu^2)*beta ) - (2/3)*log(beta);



%a = (bm/2)/sqrt(2/pi);
%Ju = (beta^2)/(2*(a^2)) - 2*log(beta);

%Ju = 0;

%     theta = 0.3;
%     a = 0;
%     b = 2*bm;
%
%     mu = bm;
%     v = Di;
%     aa = (a+(-1).*b).^(-1).*(a+b+a.*theta+(-1).*b.*theta+(-2).*mu).*((-1).*a.^2+(-4) ...
%         .*a.*b+(-1).*b.^2+(-6).*v+a.^2.*theta+(-2).*a.*b.*theta+b.^2.*theta+6.*a.*mu+ ...
%         6.*b.*mu+(-6).*mu.^2).*((-3).*a.^2+(-6).*a.*b+(-3).*b.^2+2.*a.^2.*theta+ ...
%         8.*a.*b.*theta+2.*b.^2.*theta+12.*v.*theta+a.^2.*theta.^2+(-2).*a.*b.*theta.^2+b.^2.* ...
%         theta.^2+12.*a.*mu+12.*b.*mu+(-12).*a.*theta.*mu+(-12).*b.*theta.*mu+(-12).*mu.^2+ ...
%         12.*theta.*mu.^2).^(-1);
%
%
%     bb = (a+(-1).*b).^(-1).*((-1).*a+(-1).*b+a.*theta+(-1).*b.*theta+2.*mu).*((-1).* ...
%         a.^2+(-4).*a.*b+(-1).*b.^2+(-6).*v+a.^2.*theta+(-2).*a.*b.*theta+b.^2.*theta+ ...
%         6.*a.*mu+6.*b.*mu+(-6).*mu.^2).*((-3).*a.^2+(-6).*a.*b+(-3).*b.^2+2.* ...
%         a.^2.*theta+8.*a.*b.*theta+2.*b.^2.*theta+12.*v.*theta+a.^2.*theta.^2+(-2).*a.*b.* ...
%         theta.^2+b.^2.*theta.^2+12.*a.*mu+12.*b.*mu+(-12).*a.*theta.*mu+(-12).*b.*theta.*mu+( ...
%         -12).*mu.^2+12.*theta.*mu.^2).^(-1);
%
%
%     Ju = -log( (theta*gamma(aa+bb)*((beta-a)^(aa-1))*((b-beta)^(bb-1)))/(gamma(aa)*gamma(bb)*((b-a)^(aa+bb+1))) + (1-theta)/(b-a) );
%
%     if beta > b
%         Ju = inf;
%     end


%     xm = 0.5;
%     alpha = bm/(bm-xm);
%     Jf = Jf + (alpha+1)*log(beta);

%k = bm;
%Jf = Jf - (k/2-1)*log(beta)+beta/2;

%lambda = 1/bm;
%Jf = Jf + lambda*beta;

%     Jf = Jf + (1/2)*(beta - bm).'/Di*(beta - bm);


% Beta Prime
%alphap = (bm * (Di + bm + bm^2))/Di;
%betap  = (2*Di + bm + bm^2)/Di;
%Jf = Jf + (sum(-(alphap - 1)*log(beta) + (alphap + betap)*log(1+beta)));

% Nothing
%Jf = 0;



% Xam = mean(Xa, 2);
%
% Aa = Xa - Xam*ones(1, ensN);
%
% Aa = inflation*Aa;
%
%
% Xa = Xam*ones(1, ensN) + Aa;

Jo = 0;

for asd = 1:ensN
    gt = gs(:, asd);
    
    Jot = (gt.'/R)*gt;
    
    %g = G(:, asd);
    %Jo1 = 1/2*g.'/R*g;
    
    Jo = Jo + Jot;
    %Jo = Jo + 1/2*G(:, asd).'/R*G(:, asd);
end
Jo = 1/2*(Jo);


% gf =  Y*ones(1, ensN) - Xa(H, :);
% Jo = (1/2)*norm((gf.'/R)*gf);

% G = dm - W*dm;
%
% Jo = 1/2*G.'/R*G;

%da = mean(Ys - Xa(H, :), 2);
%Jo = 1/2*norm(da.'/R*da);


%Jo = 1/2*norm(mean(G'/R*G, 2));

% Jf = sqrt(det(Pfr)) * 1/2*norm(d'*(S\W)*d);
% Jo = sqrt(det(R))   * 1/2*norm(G'/R*G);

% M = length(Pfr);
% Jf = (det(Pfr))^((M+1)/2) * 1/2*norm(d'*(S\W)*d);
% Jo = (det(R))^((M+1)/2)   * 1/2*norm(G'/R*G);

J = Jf + Jo + Ju;

end

