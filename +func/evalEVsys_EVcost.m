function EVsys_EVcost = evalEVsys_EVcost( decRules,cpmProdProbVec,Pxd,EVcost )

Ncomp = size( Pxd );
Ndec = size(decRules,1);
EVsys_EVcost = zeros( Ndec,2 );
for dd = 1:Ndec
    dec_d = decRules(dd,:);
    %
    EVsys_d = cpmProdProbVec; EVcost_d = 0;
    for nn = 1:Ncomp
        EVsys_d = exp( log(EVsys_d) + log(Pxd{nn}(:,dec_d(nn))) );
        EVcost_d = EVcost_d + EVcost{nn}( dec_d(nn) );
    end
    EVsys_EVcost(dd,:) = [sum(EVsys_d) EVcost_d];
end