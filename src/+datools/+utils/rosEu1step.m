function y = rosEu1step(f, tspan, y0, varargin)

defaultFDtype = "forward";

p = inputParser;
addParameter(p, "Jacobian", [])
addParameter(p, "FiniteDifferenceType", defaultFDtype)
parse(p, varargin{:});

h = tspan(2);
feval = f(tspan(1), y0);

if ~isempty(p.Results.Jacobian)

    n = length(y0);
    Jeval = decomposition(eye(n) - h*p.Results.Jacobian(tspan(1), y0), 'lu');
    y = y0 + h*(Jeval\feval);

else

    switch p.Results.FiniteDifferenceType
        case "forward"
            jvpfh = @(w) w - h*forwardjvp(w, f, tspan(1), y0, feval);
        case "backward"
            jvpfh = @(w) w - h*backwardjvp(w, f, tspan(1), y0, feval);
        case "central"
            jvpfh = @(w) w - h*centraljvp(w, f, tspan(1), y0);
        case "complex"
            jvpfh = @(w) w - h*complexjvp(w, f, tspan(1), y0);
    end

    [jvp, ~] = gmres(jvpfh, feval, [], 1e-3, []);
    y = y0 + h*jvp;

end

end

function Jvp = forwardjvp(w, f, t, y, feval)

h = sqrt(eps);
fplus1 = f(t, y + h*w);
Jvp = (fplus1 - feval)/h;

end

function Jvp = backwardjvp(w, f, t, y, feval)

h = sqrt(eps);
fminus1 = f(t, y - h*w);
Jvp = (feval - fminus1)/h;

end

function Jvp = centraljvp(w, f, t, y)

h = (eps/2)^(1/3);
hw = h*w;

fplus1 = f(t, y + hw);
fminus1 = f(t, y - hw);

Jvp = (fplus1 - fminus1)/(2*h);

end

function Jvp = complexjvp(w, f, t, y)

h = sqrt(eps);
Jvp = imag( f(t, y + (1i*h*w))/h );

end
