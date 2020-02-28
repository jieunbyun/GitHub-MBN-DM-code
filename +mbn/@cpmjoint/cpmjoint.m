classdef cpmjoint
% joint PMF: variables = [3 1 2] -> P(X3,X1,X2)
    
    properties
        scope;
        C;
        p;
    end
    
    methods
        function cpmj = cpmjoint(scope, C, p)
            % cpmjoint(<scope>,<C>,<p>)
            % eg p1 = cpmjoint( [3 1 2],[1 1 0;1 2 1;1 2 2;2 0 2],[0.5 0.3 0.2 1])
            % or p1 = cpmjoint; p1.scope = [3 1 2]; p1.C = [1 1 0;1 2 1;1 2 2;2 0 2]; p1.p = [0.5 0.3 0.2 1];
            if nargin > 0
                if ~isnumeric(scope)
                    error( 'Scope must be a numerical vector' )
                else
                    cpmj.scope = scope(:)';
                end
                
                if nargin>1
                    if ~isnumeric(C)
                        error('Event matrix must be a numerical array')
                    elseif size(C,2) ~= length(scope)
                        error( 'The # of cols in C must be same with the # of vars in scope' )
                    else
                        cpmj.C = C;
                    end
                    
                    if nargin >2
                        if ~isnumeric(p)
                            error( 'Prob vector must be a numerical vector' )
                        elseif ~isvector(p)
                            error('Prob vector must be a vector')
                        elseif size(C,1)~=length(p)
                            error('The # of rules in C must be same with that in p')
                        else
                            cpmj.p = p(:); 
                        end
                        
                    end
                end
            end
        end
    end
end