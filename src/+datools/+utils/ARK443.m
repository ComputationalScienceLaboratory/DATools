function [tspan, y] = ARK443(f, tspan, y0, varargin)


p = inputParser;
p.KeepUnmatched = true;
addParameter(p, 'RelTol', 1e-3);
addParameter(p, 'AbsTol', 1e-6);
%addParameter(p, 'InitialStepSize', []);
parse(p, varargin{:});
s = p.Results;

reltol = s.RelTol;
abstol = s.AbsTol;
%dt     = s.InitialStepSize;

%if isempty(dt)
f1 = f(tspan(1), y0);
dt = norm(y0) / norm(f1) * ...
    0.01;
%end


v = p.Unmatched;
p = inputParser;
addParameter(p, 'JacobianTime', []);
parse(p, v);
s = p.Results;
dfdt = s.JacobianTime;


tc = tspan(1);

y1c = y0;
y2c = dt * f1;
if isempty(dfdt)
    y3c = dt * (f(tspan(1), y0+dt*f1) - f1);
else
    y3c = (dt^2) * dfdt(tc, y1c);
end

a21 = 1 / 18;
a31 = 1 / 18;
a32 = 1;
a41 = 1 / 8;
a42 = 3 / 8;
a43 = 3 / 8;


c1 = 1;
c2 = 1 / 3;
c3 = 2 / 3;
c4 = 1;

beta1 = -1 / 2;
beta2 = 3 / 2;
beta3 = -3 / 2;
beta4 = 2;
beta0 = -3 / 2;

u12 = 1;
u13 = 1 / 2;
u22 = 5 / 18;
u23 = 0;
u32 = -7 / 18;
u33 = -1 / 6;
u42 = 1 / 8;
u43 = 0;

b1h = 0.194;
%b1h = -1/4;
b2h = 3 * b1h;
b3h = 3 / 4 - b2h;
b0h = 1 - b1h - b2h - b3h;

y = zeros(numel(y0), numel(tspan));

y(:, 1) = y0;

for sl = 2:numel(tspan)

    y1c = y(:, sl-1);

    tc = tspan(sl - 1);
    tend = tspan(sl);

    if tc + dt > tend
        dtnew = tend - tc;

        y2c = (dtnew / dt) * y2c;
        y3c = ((dtnew / dt)^2) * y3c;

        dt = dtnew;
    end

    while tc < tend

        %dt

        Y1 = y1c + u12 * y2c + u13 * y3c;
        fY1 = f(tc+c1*dt, Y1);
        Y2 = y1c + (dt * a21) * fY1 + u22 * y2c + u23 * y3c;
        fY2 = f(tc+c2*dt, Y2);
        Y3 = y1c + (dt * a31) * fY1 + (dt * a32) * fY2 + u32 * y2c + u33 * y3c;
        fY3 = f(tc+c3*dt, Y3);
        Y4 = y1c + (dt * a41) * fY1 + (dt * a42) * fY2 + (dt * a43) * fY3 + u42 * y2c + u43 * y3c;


        y1cn = Y4;
        y1hat = y1c + (dt * b1h) * fY1 + (dt * b2h) * fY2 + (dt * b3h) * fY3 + b0h * y2c;

        sc = abstol + max(abs(y1cn), abs(y1hat)) * reltol;

        err = rms((y1cn - y1hat)./sc);

        orderE = 3;
        fac = 0.38^(1 / (orderE + 1));
        %fac = 0.9;


        facmin = 0.2;
        if err > 1 || isnan(err)
            % Reject timestep
            facmax = 1;

            % adjust step-size
            dtnew = dt * min(facmax, max(facmin, fac*(1 / err)^(1 / (orderE + 1))));
        else

            % Accept time step
            tc = tc + dt;
            y1c = y1cn;
            fY4 = f(tc+c4*dt, Y4);
            %[gY4, fY4] = dfdt(tc + c4*dt, Y4);
            %y3c = (dt^2)*gY4;
            y3c = dt * beta1 * fY1 + dt * beta2 * fY2 + dt * beta3 * fY3 + dt * beta4 * fY4 + beta0 * y2c;
            y2c = dt * fY4;

            facmax = 3;

            % adjust step-size
            dtnew = dt * min(facmax, max(facmin, fac*(1 / err)^(1 / (orderE + 1))));

            if tc + dtnew > tend && tc < tend
                dtnew = tend - tc;
            end

        end

        y2c = (dtnew / dt) * y2c;
        y3c = ((dtnew / dt)^2) * y3c;

        dt = dtnew;

        %y3c = (dt^2)*dfdt(tc, y1c);

    end

    y(:, sl) = y1c;

end

y = y.';

end
