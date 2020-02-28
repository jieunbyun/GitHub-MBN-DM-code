function EVcost_basis = evalEVcost( basisDecRule,EVcost )

Ncomp = length(basisDecRule);

EVcost_basis = 0;
for nn = 1:Ncomp
    EVcost_basis = EVcost_basis + EVcost{ nn }( basisDecRule(nn) );
end