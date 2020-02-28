function proxEVsys = getProxEVSys( Cpm,var,basisDecRule,Pxd,CpmProd )

import mbn.*
import func.*

Ncomp = length( var.X );
Ndec = cellfun( @(x) length( x.C ), Cpm(var.D) );

EVsys_basis = CpmProd.p;
for nn = 1:Ncomp
    EVsys_basis = exp( log(EVsys_basis) + log(Pxd{nn}(:,basisDecRule(nn))) );
end

proxEVsys = {};
for nn = 1:Ncomp
    proxEVsys_n = zeros( Ndec(nn),1 );
    for dd = 1:Ndec(nn)
        if dd == basisDecRule(nn)
            proxEVsys_n(dd) = sum( EVsys_basis );
        else
            EVsys_basis_nd = exp( log(EVsys_basis) - log( Pxd{nn}(:,basisDecRule(nn))) + log(Pxd{nn}(:,dd)) );
            proxEVsys_n(dd) = sum( EVsys_basis_nd );
        end
    end
    proxEVsys = [proxEVsys; {proxEVsys_n}];
end

