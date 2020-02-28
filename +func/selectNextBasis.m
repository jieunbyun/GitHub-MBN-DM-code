function [basisDecRule,nextBasisId] = selectNextBasis( decOpt_EVsys,decOpt_prox_EVprox,decOpt_prox )

[decOpt_EVsys_sort, decOpt_EVsys_sortId] = sort( decOpt_EVsys );
decOpt_EVsys_prox_sort = decOpt_prox_EVprox( decOpt_EVsys_sortId );
div = -( decOpt_EVsys_prox_sort(2:end) - decOpt_EVsys_prox_sort(1:end-1) + eps )./( decOpt_EVsys_sort(2:end) - decOpt_EVsys_sort(1:end-1) +eps );
[~,nextBasisId_sort] = max( div );
nextBasisId = decOpt_EVsys_sortId( nextBasisId_sort+1 );
basisDecRule = decOpt_prox( nextBasisId,: );
       