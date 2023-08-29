function [tspan, y, l] = dp54tlm(f, tspan, y0, l0, options)

%% Get options
if nargin < 4
    options = struct;
end

reltol = odeget(options, 'RelTol', 1e-3);
abstol = odeget(options, 'AbsTol', 1e-6);
J = odeget(options, 'Jacobian', ...
    @(t, y) otp.utils.derivatives.jacobian(f, t, y, 'FD'));

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
if N > 1
    error('Unsupported');
end
tc = tspan(1)*ones(1, N);
yc = y0;
lc = l0;

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

y = yc;

K7 = J(tc, yc);

hold = h;

for tcheck = 2:numel(tspan)

    tc = tspan(tcheck - 1);
    tend = tspan(tcheck);

    %% Initialize time step
    h = hold;

    %% Initialize finalizing booleans
    endit = false;
    ended = false;

    %% Perform the time stepping
    while ~ended

        if ((tc + h) >= tend)
            % save the timestep
            hold = h;
            % set the ending time step to something that will hit our target
            h = tend - tc;
            endit = true;
        end
        
        hK1 = h.*K7;
        hk1 = h.*k7;

        y2 = yc + a21*hk1;
        hk2 = h.*f(tc + c2*h, yc + a21*hk1);

        y3 = yc + a31*hk1 + a32*hk2;
        hk3 = h.*f(tc + c3*h, y3);

        y4 = yc + a41*hk1 + a42*hk2 + a43*hk3;
        hk4 = h.*f(tc + c4*h, y4);

        y5 = yc + a51*hk1 + a52*hk2 + a53*hk3 + a54*hk4;
        hk5 = h.*f(tc + c5*h, y5);

        y6 = yc + a61*hk1 + a62*hk2 + a63*hk3 + a64*hk4 + a65*hk5;
        hk6 = h.*f(tc + h, y6);

        y1cn = yc + a71*hk1 + a73*hk3 + a74*hk4 + a75*hk5 + a76*hk6;

        hk7 = h.*f(tc + h, y1cn);
        y1hat = yc + bb1*hk1 + bb3*hk3 + bb4*hk4 + bb5*hk5 + bb6*hk6 + bb7*hk7;

        sc = abstol + max(abs(yc), abs(y1cn))*reltol;
        err = sqrt(mean( ((y1cn - y1hat)./sc).^2, 1 ));

        % reject
        if (err > 1) && (~ended)
            h = h.*min(0.9, max(facmin, fac*(1./err).^(1/5)));
            % don't end that ensemble run if we rejected it
            endit = false;
            continue;
        end

        k7 = hk7./h;
        tc= tc + h;
        yc = y1cn;

        % compute TLM
        l2 = lc + a21*hK1;
        hK2 = h.*J(tc + c2*h, y2)*l2;
        l3 = lc + a31*hK1 + a32*hK2;
        hK3 = h.*J(tc + c3*h, y3)*l3;
        l4 = lc + a41*hK1 + a42*hK2 + a43*hK3;
        hK4 = h.*J(tc + c4*h, y4)*l4;
        l5 = lc + a51*hK1 + a52*hK2 + a53*hK3 + a54*hK4;
        hK5 = h.*J(tc + c5*h, y5)*l5;
        l6 = lc + a61*hK1 + a62*hK2 + a63*hK3 + a64*hK4 + a65*hK5;
        hK6 = h.*J(tc + h, y6)*l6;
        l7 = lc + a71*hK1 + a73*hK3 + a74*hK4 + a75*hK5 + a76*hK6;
        K7 = J(tc + h, y1cn)*l7;

        h = h...
            .*min(facmax, max(facmin, fac*(1./err).^(1/5)));

        % end what we need to end and set the time step to zero
        if endit || isnan(h)
            ended = true;
            h = 0;
        end

    end

    y(:, :, tcheck - 1) = yc;

end

l = lc;
