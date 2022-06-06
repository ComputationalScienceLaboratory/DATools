function [tspan, y] = RODAS4(f, tspan, y0, varargin)

% coeffA = zeros(4);
% coeffG = zeros(4);

gamma = 1 / 4;

a1 = 0.0;
a2 = 0.386;
a3 = 0.210;
a4 = 0.630;
a5 = 1.0;
a6 = 1.0;

c21 = -0.5668800000000000e+1;
c31 = -0.2430093356833875e+1;
c32 = -0.2063599157091915e+0;
c41 = -0.1073529058151375e+0;
c42 = -0.9594562251023355e+1;
c43 = -0.2047028614809616e+2;
c51 = 0.7496443313967647e+1;
c52 = -0.1024680431464352e+2;
c53 = -0.3399990352819905e+2;
c54 = 0.1170890893206160e+2;
c61 = 0.8083246795921522e+1;
c62 = -0.7981132988064893e+1;
c63 = -0.3152159432874371e+2;
c64 = 0.1631930543123136e+2;
c65 = -0.6058818238834054e+1;

a21 = 0.1544000000000000e+1;
a31 = 0.9466785280815826e+0;
a32 = 0.2557011698983284e+0;
a41 = 0.3314825187068521e+1;
a42 = 0.2896124015972201e+1;
a43 = 0.9986419139977817e+0;
a51 = 0.1221224509226641e+1;
a52 = 0.6019134481288629e+1;
a53 = 0.1253708332932087e+2;
a54 = -0.6878860361058950e+0;
a61 = a51;
a62 = a52;
a63 = a53;
a64 = a54;
a65 = 1.0;

m1 = a51;
m2 = a52;
m3 = a53;
m4 = a54;
m5 = 1;
m6 = 1;

mh1 = a51;
mh2 = a52;
mh3 = a53;
mh4 = a54;
mh5 = 1;

p = inputParser;
addParameter(p, 'RelTol', 1e-3);
addParameter(p, 'AbsTol', 1e-6);
jvpe = sqrt(eps);
addParameter(p, 'JacobianVectorProduct', @(t, y, u) (f(t, y+jvpe*u) - f(t, y))/jvpe);
addParameter(p, 'OutputFcn', []);
addParameter(p, 'InitialStepSize', 1e-3);

parse(p, varargin{:});
s = p.Results;

reltol = s.RelTol;
abstol = s.AbsTol;
Jvpa = s.JacobianVectorProduct;
df = s.OutputFcn;
dt = s.InitialStepSize;


yc = y0;


tc = tspan(1);
tend = tspan(end);

if tc + dt > tend
    dt = tend - tc;
end

restarts = [];

orderE = 3;

lnsoltol = reltol;

s = 1;

while tc < tend

    dg = dt * gamma;

    A = @(x) x / dg - Jvpa(tc, yc, x);

    t1 = tc + dt * a1;
    Y1 = yc;
    fe1 = f(t1, Y1);
    u1rhs = fe1;

    [u1, ~] = gmres(A, u1rhs, restarts, lnsoltol);


    t2 = tc + dt * a2;
    Y2 = yc + a21 * u1;
    fe2 = f(t2, Y2);


    u2rhs = ((fe2 + c21 * u1 / dt));

    [u2, ~] = gmres(A, u2rhs, restarts, lnsoltol);


    t3 = tc + dt * a3;
    Y3 = yc + a31 * u1 + a32 * u2;
    fe3 = f(t3, Y3);
    u3rhs = ((fe3 + (c31 * u1 + c32 * u2) / dt));

    [u3, ~] = gmres(A, u3rhs, restarts, lnsoltol);


    t4 = tc + dt * a4;
    Y4 = yc + a41 * u1 + a42 * u2 + a43 * u3;
    fe4 = f(t4, Y4);
    u4rhs = ((fe4 + (c41 * u1 + c42 * u2 + c43 * u3) / dt));

    [u4, ~] = gmres(A, u4rhs, restarts, lnsoltol);


    t5 = tc + dt * a5;
    Y5 = yc + a51 * u1 + a52 * u2 + a53 * u3 + a54 * u4;
    fe5 = f(t5, Y5);
    u5rhs = ((fe5 + (c51 * u1 + c52 * u2 + c53 * u3 + c54 * u4) / dt));

    [u5, ~] = gmres(A, u5rhs, restarts, lnsoltol);


    t6 = tc + dt * a6;
    Y6 = yc + a61 * u1 + a62 * u2 + a63 * u3 + a64 * u4 + a65 * u5;
    fe6 = f(t6, Y6);
    u6rhs = ((fe6 + (c61 * u1 + c62 * u2 + c63 * u3 + c64 * u4 + c65 * u5) / dt));

    [u6, ~] = gmres(A, u6rhs, restarts, lnsoltol);


    yhat = yc + mh1 * u1 + mh2 * u2 + mh3 * u3 + mh4 * u4 + mh5 * u5;
    ycn = yc + m1 * u1 + m2 * u2 + m3 * u3 + m4 * u4 + m5 * u5 + m6 * u6;


    sc = abstol + max(abs(ycn), abs(yhat)) * reltol;

    err = rms((ycn - yhat)./sc);

    %fac = 0.38^(1/(orderE + 1));
    fac = 0.9;

    facmin = 0.2;
    if err > 1 || isnan(err)
        % Reject timestep
        facmax = 1;

        % adjust step-size
        dt = dt * min(facmax, max(facmin, fac*(1 / err)^(1 / (orderE + 1))));
    else
        % Accept time step
        tc = tc + dt;
        yc = ycn;


        facmax = 6;

        % adjust step-size
        dt = dt * min(facmax, max(facmin, fac*(1 / err)^(1 / (orderE + 1))));

        if tc + dt > tend
            dt = tend - tc;
        end

        %t(s + 1) = tc;
        %y(:, s + 1) = yc;

        s = s + 1;

        if ~isempty(df)
            df(tc, yc, []);
        end


    end


end

y = [y0, yc].';

end
