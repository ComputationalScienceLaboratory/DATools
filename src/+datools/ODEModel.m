classdef ODEModel < handle
%ODEModel This class handles ODEmodel
%   The user can use OTP or provide relevant information   

    properties
        F           % RHS of the ODE Model
        J           % Jacobian
        JVP         % Jacobian Vector Product
        JAVP        % Jacobian Adjoint Vector Product
        HVP         % Hessian Vector Product
        HVAP        % Hessian Vector Adjoint Product
        TimeSpan    % Time Span for One Step Evolution
        Y0          % Initial state of the system
        NumVars     % Number of State Variables
    end
    
    methods
        function obj = ODEModel(varargin)
           p = inputParser;
           p.KeepUnmatched = true;
           addParameter(p, 'OTPObject', [])
           addParameter(p, 'F', @(x) x);
           addParameter (p, 'J', []);
           addParameter(p, 'JVP', []);
           addParameter(p, 'JAVP', []);
           addParameter(p, 'HVP', []);
           addParameter(p, 'HVAP', []);
           addParameter(p, 'TimeSpan', []);
           addParameter(p, 'Y0', []);
           addParameter(p, 'NumVars', []);
           
           parse(p, varargin{:});
           
           s = p.Results;
           
           OTPObject = s.OTPObject;
           if isa(OTPObject, 'otp.Problem')
               % cal the desired function
               obj.assignValuesOTP(OTPObject);
           else
               obj.F = s.F;
               obj.J = s.J;
               obj.JVP = s.JVP;
               obj.JAVP = s.JAVP;
               obj.HVP= s.HVP;
               obj.HVAP = s.HVAP;
               obj.TimeSpan = s.TimeSpan;
               obj.Y0 = s.Y0;
               obj.NumVars = s.NumVars;
           end
        end
        
        
        function assignValuesOTP(obj, ode)
            % put a check here
            obj.F = ode.RHS.F;
            obj.J = ode.RHS.Jacobian;
            obj.JVP = ode.RHS.JacobianVectorProduct;
            obj.JAVP = ode.RHS.JacobianAdjointVectorProduct;
            obj.HVP = ode.RHS.HessianVectorProduct;
            obj.HVAP = ode.RHS.HessianAdjointVectorProduct;
            obj.TimeSpan = ode.TimeSpan;
            obj.Y0 = ode.Y0;
            obj.NumVars = ode.NumVars;
        end
        
        function assignValuesUser()
            % complete later
            % provide sufficient checks
        end
        
    end
end