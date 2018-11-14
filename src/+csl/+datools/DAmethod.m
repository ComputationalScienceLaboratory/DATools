classdef DAmethod < handle
    %DAMETHOD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Problem
        Solver
        H                % H is our observation operator
        CurrentBestGuess % This is our current best guess of the system
    end
    
    methods
        function obj = DAmethod(problem, solver, H)
            obj.Problem = problem;
            obj.Solver  = solver;
            if (nargin > 2)
                obj.H   = H;
            else
                obj.H   = 1;
            end
            obj.CurrentBestGuess = obj.Problem.Y0;
        end
        function forecast(obj)
            % This is the basic implementation of the forecast. I would expect
            % most data assimilation methods to modify this (for instance
            % ensemble methods will propogate an ensemble.
            
            % Use our ode solver to solve the IVP
            [t, y] = obj.Solver(obj.Problem.F, obj.Problem.TimeSpan, obj.Problem.Y0);
            
            y = y'; % We are all sane people here.
            
            % Set the new timespan to be that of the end time + dt
            deltat = diff(obj.Problem.TimeSpan);
            obj.Problem.TimeSpan = [t(end), t(end)+deltat];
            
            % Set Y0 to the last iteration of the forecast
            obj.Problem.Y0 = y(:,end);
            obj.CurrentBestGuess = y(:, end);
        end
    end
    
end

