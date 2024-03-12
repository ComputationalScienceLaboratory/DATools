classdef Polynomial
    

    properties
        p
    end

    methods
        function obj = Polynomial(p)
            obj.p = p;
        end

        function obj = minus(obj, obj2)
            if ~isa(obj,"datools.utils.Polynomial")
                obj = datools.utils.Polynomial(obj);
            end

            if ~isa(obj2,"datools.utils.Polynomial")
                obj2 = datools.utils.Polynomial(obj2);
            end
            p1 = obj.p;
            p2 = obj2.p;

            n1 = numel(p1);
            n2 = numel(p2);
            p1 = padarray(p1,[0, max(n2 - n1, 0)],0, "pre");
            p2 = padarray(p2,[0, max(n1 - n2, 0)],0, "pre");
            
            pnew = p1 - p2;
            obj.p = pnew;
        end

        function obj = mtimes(obj, obj2)
            if ~isa(obj,"datools.utils.Polynomial")
                obj = datools.utils.Polynomial(obj);
            end
            if ~isa(obj2,"datools.utils.Polynomial")
                obj2 = datools.utils.Polynomial(obj2);
            end
            p1 = obj.p;
            p2 = obj2.p;
            pnew = conv(p1, p2);
            obj.p = pnew;
        end
    end

end
