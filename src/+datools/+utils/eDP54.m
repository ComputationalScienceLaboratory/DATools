function [tspan, y] = eDP54(f, tspan, y0, options)

%% Get options
if nargin < 4
    options = struct;
end

reltol = odeget(options, 'RelTol', 1e-3);
abstol = odeget(options, 'AbsTol', 1e-6);

%% Define error control factors
facmax = 6;
facmin = 0.2;
fac = 0.38^(1/(4 + 1));

%% Butcher Tableau
a21 = 1/5;
a31 = 3/40;
a32 = 9/40;
a41 = 44/45;
a42 = -56/15;
a43 = 32/9;
a51 = 19372/6561;
a52 = -25360/2187;
a53 = 64448/6561;
a54 = -212/729;
a61 = 9017/3168;
a62 = -355/33;
a63 = 46732/5247;
a64 = 49/176;
a65 = -5103/18656;
a71 = 35/384;
% a72 = 0;
a73 = 500/1113;
a74 = 125/192;
a75 = -2187/6784;
a76 = 11/84;

bb1 = 5179/57600;
% bb2 = 0;
bb3 = 7571/16695;
bb4 = 393/640;
bb5 = -92097/339200;
bb6 = 187/2100;
bb7 = 1/40;

% c1 = 0;
c2 = 1/5;
c3 = 3/10;
c4 = 4/5;
c5 = 8/9;
% c6 = 1;
% c7 = 1;

%% Define time and state
[n, N] = size(y0);
tc = tspan(1)*ones(1, N);
yc = y0;

%% Initial step size selection
k7 = f(tc, yc);
sc = abstol + abs(y0)*reltol;
d0 = sqrt(mean( (y0./sc).^2, 1 ));
d1 = sqrt(mean( (k7./sc).^2, 1 ));
h0 = 0.01*(d0./d1);
y1 = y0 + h0.*k7;
fy1 = f(tc + h0, y1);
d2 = sqrt(mean( ((k7 - fy1)./sc).^2, 1 ))./h0;
h1 = (0.01./max(d1, d2)).^(1/6);
h = min(100*h0, h1);

y = zeros(n, N, numel(tspan) - 1, 'like', k7);

hold = h;

for tcheck = 2:numel(tspan)

    tc = tspan(tcheck - 1)*ones(1, N);
    tend = tspan(tcheck);

    %% Initialize time step
    h = hold;

    %% Initialize finalizing booleans
    endit = false(1, N);
    ended = false(1, N);

    %% Perform the time stepping
    while any(~ended)

        indsendcandidate = (~ended) & ((tc + h) >= tend);
        % save the timestep
        hold(indsendcandidate) = h(indsendcandidate);
        % set the ending time step to something that will hit our target
        h(indsendcandidate) = tend - tc(indsendcandidate);
        endit(indsendcandidate) = true;

        hk1 = h.*k7;
        hk2 = h.*f(tc + c2*h, yc + a21*hk1);
        hk3 = h.*f(tc + c3*h, yc + a31*hk1 + a32*hk2);
        hk4 = h.*f(tc + c4*h, yc + a41*hk1 + a42*hk2 + a43*hk3);
        hk5 = h.*f(tc + c5*h, yc + a51*hk1 + a52*hk2 + a53*hk3 + a54*hk4);
        hk6 = h.*f(tc + h, yc + a61*hk1 + a62*hk2 + a63*hk3 + a64*hk4 + a65*hk5);

        y1cn = yc + a71*hk1 + a73*hk3 + a74*hk4 + a75*hk5 + a76*hk6;

        hk7 = h.*f(tc + h, y1cn);
        y1hat = yc + bb1*hk1 + bb3*hk3 + bb4*hk4 + bb5*hk5 + bb6*hk6 + bb7*hk7;

        sc = abstol + max(abs(yc), abs(y1cn))*reltol;
        err = sqrt(mean( ((y1cn - y1hat)./sc).^2, 1 ));

        % reject
        indsreject = (err > 1) & (~ended);
        h(indsreject) = h(indsreject)...
            .*min(0.9, max(facmin, fac*(1./err(indsreject)).^(1/5)));
        % don't end that ensemble run if we rejected it
        endit(indsreject) = false;

        % accept
        indsaccept = (err <= 1) & (~ended);

        % I have to do it this way because matlab indexing is bad sometimes
        k7new = hk7./h;
        k7(:, indsaccept) = k7new(:, indsaccept);
        tc(indsaccept) = tc(indsaccept) + h(indsaccept);
        yc(:, indsaccept) = y1cn(:, indsaccept);

        h(indsaccept) = h(indsaccept)...
            .*min(facmax, max(facmin, fac*(1./err(indsaccept)).^(1/5)));

        % end what we need to end and set the time step to zero
        ended(endit) = true;
        h(ended) = 0;

        ended(isnan(h)) = true;

    end

    y(:, :, tcheck - 1) = yc;

end
