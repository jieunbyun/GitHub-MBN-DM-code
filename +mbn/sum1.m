function cpmnew = sum1( cpm,X,Xstay )
%{
Sum operation of a cpm (cpm) over given variables (X)
Input:
    cpm: a cpm
    X: Nx x 1 array of variable
    <Xstay>: a scalar of sum option (0-default: X being summed up; 1: other than X being summed up)
Output:
    cpmnew: a cpm after sum operation

Ex: 
cpm = cpmcond([3 5 1 2],2,[1 1 1 3;2 2 1 3;1 1 2 1;1 2 2 1;2 1 2 2],[.4 .6 .3 .7 1]');
cpm = cpmjoint([3 5 1 2],[1 1 1 3;2 2 1 3;1 1 2 1;1 2 2 1;2 1 2 2],[.2 .1 .25 .4 .05]');

X = [5 2 4];
cpm2 = sum1( cpm,X );
cpm2 = sum1( cpm,X,1 );
%}

import mbn.*
if isa(cpm,'cpmcond')
    X = setdiff(X,cpm.scope((cpm.nc+1):end),'stable'); % conditioned var's are not summed up.
end

if nargin>2
   if Xstay
       [Xin,Xi_idc] = intersect(cpm.scope,X);
       if isa(cpm,'cpmcond')
           Xin = [Xin(:)' cpm.scope( (cpm.nc+1):end )];
           Xi_idc = [Xi_idc(:)' (cpm.nc+1):length(cpm.scope)];
       end       
   else
       [Xin,Xi_idc] = setdiff(cpm.scope,X);
   end
else
    [Xin,Xi_idc] = setdiff(cpm.scope,X);
end

Cnew = []; pnew = [];
Cold = cpm.C(:,Xi_idc); pold = cpm.p;
while ~isempty(pold)
   c = Cold(1,:);
   id = ismember( Cold,c,'rows' );
   Cnew = [Cnew;c];
   pnew = [pnew;sum(pold(id))];
   Cold(id,:) = [];
   pold(id) = [];
end

[~,Xi_idc2] = sort(Xi_idc);
Cnew = Cnew(:,Xi_idc2); Xin = Xin(Xi_idc2);
if isa(cpm,'cpmcond')
    cpmnew = cpmcond;
    cpmnew.nc = sum(Xi_idc<=cpm.nc);
else
    cpmnew = cpmjoint;
end
cpmnew.scope = Xin;
cpmnew.C = Cnew;
cpmnew.p = pnew;