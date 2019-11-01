function [t, y] = solverprojection(f, ts, y0, solver, projection)

[t, y] = solver(f, ts, y0);

y = (projection*(y.')).';

end