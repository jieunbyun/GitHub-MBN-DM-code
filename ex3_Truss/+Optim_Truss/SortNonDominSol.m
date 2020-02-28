function [Dopt,EVs,EC] = SortNonDominSol( D,EVs,EC )
%{
Sort non-dominated solutions (Pareto Solutions) in D

INPUT:
    D: a set of decision rules - Nsol x Nx matrix
    EVs: E[Vs|d] value - Nsol x 1 vector
    EC: cost E[C|d] (exact) - Nsol x 1 vector
OUTOUT:
    Dopt: non-dominated solutions - Nopt x Nx matrix
    EVs: E[Vs|d] for non-dominated solutions - Nopt x 1 vector
    EC: E[C|d] for non-dominated solutions - Nopt x 1 vector
%}

D2 = D; EVs2=EVs; EC2 = EC;
Dopt = []; EVs=[]; EC = [];
while size(D2,1)
    d_ = D2(1,:); evs_ = EVs2(1); ec_ = EC2(1);
    if ~all( EVs2(2:end)>evs_ | EC2(2:end)>ec_ )
        D2(1,:) = []; EVs2(1)=[]; EC2(1)=[];
    else
        Dopt=[Dopt;d_]; EVs=[EVs;evs_]; EC=[EC;ec_];
        tmp_ = find( ~(EVs2(2:end)<evs_ | EC2(2:end)<ec_) ) + 1;
        D2( [1 tmp_(:)'],: ) = []; EVs2( [1 tmp_(:)'] ) = []; EC2( [1 tmp_(:)'] ) = [];
    end
end