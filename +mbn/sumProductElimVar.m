function Cpmnew = sumProductElimVar( Cpm,elimvar,B )
%{
Eliminate a rv from a set of CPMs by sum-product

Input:
    Cpm: |M| x 1 cell array of CPM
    elimvar: scalar of random variable
Output:
    Cpmnew: |M'| x 1 cell array of CPM
%}

import mbn.*

Nm = length(Cpm);
elim = zeros( Nm,1 );
for ii = 1:Nm
    if ismember( elimvar,Cpm{ii}.scope )
        elim(ii) = 1;
    end
end

Cpmnew = Cpm(~elim);
Cpm1 = Cpm(logical(elim));
if ~isempty(Cpm1)
    cpm_ = Cpm1{1};
    for ii = 2:sum(elim)
        cpm_ = prod1( cpm_,Cpm1{ii},B );
    end
    cpm_ = sum1( cpm_,elimvar );
    Cpmnew = [Cpmnew;{cpm_}];
end