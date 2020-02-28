classdef cpmcond
% Conditional pmf: P(scopec|scopep)
    properties
        scope;
        nc; % # of scopec
        C;
        p;
    end
    
    methods
        function cpmc = cpmcond(scope, nc, C, p)
            % cpmcond(<scope>, <scopep>, <C>, <p>)
            % eg p1 = cpmcond([3 1 2],1,[1 1 0;1 2 1;1 2 2;2 0 2],[0.5 0.3 0.2 1])
            % or p1=cpmcond; p1.scope = [3 1 2]; p1.nc = 1; p1.C=[1 1 0;1 2 1;1 2 2;2 0 2]; p1.p=[0.5 0.3 0.2 1];
            if nargin>0
                if ~isnumeric(scope)
                    error( 'Scope must be a numerical vector' )
                else
                    cpmc.scope = scope(:)';
                end
                
                if nargin>1
                    if ~isscalar(nc) || nc <= 0
                        error( 'Scope2 must be a positive scalar' )
                    else
                        cpmc.nc = nc(:)';
                    end
                    
                    if nargin>2
                        if ~isnumeric(C)
                            error('Event matrix must be a numerical array')
                    elseif size(C,2) ~= length(scope)
                        error( 'The # of cols in C must be same with the # of vars in scope' )
                    else
                        cpmc.C = C;
                        end
                    
                        if nargin>3
                            if ~isnumeric(p)
                                error( 'Prob vector must be a numerical vector' )
                            elseif ~isvector(p) && ~isempty(p)
                                error('Prob vector must be a vector')
                            elseif size(C,1)~=length(p) && ~isempty(p) && ~isempty(C)
                                error('The # of rules in C must be same with that in p')
                            else
                                cpmc.p = p(:); 
                            end
                        end
                    end
                end
            end
        end
    end
end