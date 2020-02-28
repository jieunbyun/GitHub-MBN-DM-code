function [Pxd,CpmProd] = evalPxd( Cpm,var,B,varProd )

import mbn.*

Ncomp = length( var.X );
Ndec = cellfun( @(x) length( x.C ), Cpm(var.D) );

CpmProd = Cpm{varProd(1)};
for vv = 2:length(varProd)
    CpmProd = prod1( CpmProd,Cpm{varProd(vv)},B );
end
Nc_CpmProd = size( CpmProd.C,1 );

Pxd = {};
for nn = 1:Ncomp
    pX_n = zeros( Nc_CpmProd,Ndec(1) );
    Cpm_Xn = Cpm{var.X(nn)};
    comVar = intersect( CpmProd.scope,Cpm_Xn.scope );
    [~,comVarID_inCpmProd] = ismember( comVar,CpmProd.scope );
    [~,comVarID_inCpmXn] = ismember( comVar,Cpm_Xn.scope );
    [~,decID_inCpmXn] = ismember( var.D(nn),Cpm_Xn.scope );

    for cc = 1:Nc_CpmProd
        assign_c = CpmProd.C(cc,comVarID_inCpmProd);
        for dd = 1:Ndec(nn)
            compcheck_ncd = compats( Cpm_Xn.C(:,[comVarID_inCpmXn decID_inCpmXn]),[assign_c dd],B([comVar var.D(nn)]) );
            pX_n( cc,dd ) = sum( Cpm_Xn.p(compcheck_ncd) );
        end
    end
    
    Pxd = [Pxd; {pX_n}];
end

