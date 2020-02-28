function [cpmnew,B] = prod1(cpm1,cpm2,B)
%{
Product operation of two CPMs (cpm1,cpm2)
Input:
    cpm1,cpm2: a cpm 
    B: Nx x 1 cell array of basic state indicator
Output:
    cpmnew: a cpm after product operation
    B: Nx x 1 cell array of basic state indicator after product operation

Ex:
cpm1 = cpmcond([1 3 5 2],2,[1 1 1 3;2 2 1 3;2 1 2 1;2 2 2 1;1 2 2 2],[.4 .6 .1 .9 1]');
cpm2 = cpmcond([4 2 3],2,[1 1 1;2 1 1;2 2 2],[.2 .8 1]');
B{1} = eye(2); B{2,1} = [eye(2); 1 1]; B{3} = eye(2); B{4} = eye(2); B{5} = eye(2); 
[cpmnew,B] = prod1(cpm1,cpm2,B);
%}

import mbn.*

if isa(cpm1,'cpmcond')
    scope3 = cpm1.scope(1:cpm1.nc);
    if isa(cpm2,'cpmcond')
        if ~isempty( intersect(scope3,cpm2.scope(1:cpm2.nc)) )
            error('Child variables must not be overlapped')
        else
            scope3 = [scope3 cpm2.scope(1:cpm2.nc)];
            scope3p = [cpm1.scope((cpm1.nc+1):end) cpm2.scope((cpm2.nc+1):end)];
        end
    else
        if ~isempty( intersect(scope3,cpm2.scope) )
            error('Child variables must not be overlapped')
        else
            scope3 = [scope3 cpm2.scope];
            scope3p = cpm1.scope((cpm1.nc+1):end);            
        end
    end
else
    scope3 = cpm1.scope;
    if isa(cpm2,'cpmcond')
        if ~isempty( intersect(scope3,cpm2.scope(1:cpm2.nc)) )
            error('Child variables must not be overlapped')
        else
            scope3 = [scope3 cpm2.scope(1:cpm2.nc)];
            scope3p = cpm2.scope((cpm2.nc+1):end);
        end
    else
        if ~isempty( intersect(scope3,cpm2.scope) )
            error('Child variables must not be overlapped')
        else
            scope3 = [scope3 cpm2.scope]; 
            scope3p = [];
        end
    end
end
                        
if isempty( setdiff(scope3p,scope3) )
    cpmnew = cpmjoint;
else
    cpmnew = cpmcond;
    cpmnew.nc = length( unique(scope3) );
    scope3 = [scope3 scope3p];
end
scope3 = unique(scope3,'stable');

[Xp,Xpid1,Xpid2] = intersect(cpm1.scope,cpm2.scope);
[~,Xpid3] = intersect( scope3,Xp );
[X1,X1id1] = setdiff(cpm1.scope,Xp);
[~,X1id3] = intersect( scope3,X1 );
[X2,X2id2] = setdiff(cpm2.scope,Xp);
[~,X2id3] = intersect( scope3,X2 );

C3 = []; p3 = [];
Nr1 = size(cpm1.C,1); Nr2 = size(cpm2.C,1); Nx3 = length( scope3 ); Nxp = length(Xp);
for ii = 1:Nr2
    Cnew = zeros(Nr1,Nx3); Cid = 1:Nr1; c2i = cpm2.C(ii,:);
    for jj = 1:Nxp
        X1j = Xpid1(jj); X2j = Xpid2(jj); X3j = Xpid3(jj); Xj = Xp(jj);
        [ccj,cj,Bj] = compat1( cpm1.C(Cid,X1j),c2i(X2j),B{Xj} );
        
        Cid = Cid(ccj); Cnew = Cnew(ccj,:); Cnew(:,X3j) = cj; B{Xj} = Bj; 
    end
    Cnew(:,X1id3) = cpm1.C(Cid,X1id1);
    Cnew(:,X2id3) = repmat(c2i(X2id2),length(Cid),1);
    pnew = cpm1.p(Cid)*cpm2.p(ii);
    
    C3 = [C3; Cnew]; p3 = [p3;pnew];
end

cpmnew.scope = scope3; 
cpmnew.C = C3;
cpmnew.p = p3;