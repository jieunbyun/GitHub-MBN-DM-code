function cpmnew = sumProductVE( Cpm,elimvars,B )
%{
sum-product-variable-elimination

Input:
    Cpm: |M| x 1 cell array of CPM
    elimvar: scalar of random variable
Output:
    cpmnew: a cpm

%}

import mbn.*

Cpmnew = Cpm;
for ii = 1:length(elimvars)
    Cpmnew = sumProductElimVar( Cpmnew,elimvars(ii),B );
end
cpmnew = Cpmnew{1};
for ii = 2:length(Cpmnew)
    cpmnew = prod1( cpmnew,Cpmnew{ii},B );
end