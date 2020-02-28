function D = GenMut( Nmu,Nd,D,rmu )
%{
Generate Nmu decision rules by mutation from the ones in D
INPUT: 
    Nmu: # of decision rules to be generated
    Nd: # of decision alternatives
    D: a set of rules from which the new generation is produced
    rmu: ratio of genetic changes
OUTPUT:
    D: Nmu x Nx matrix with generated rules
%}

Nsol = size(D,1); Nx = size(D,2); D2 = zeros(Nmu,Nx);
for ii = 1:Nmu
    tmp_ = randsample(Nsol,1);
    tmp2_ = randsample(2,Nx,'true',[1-rmu rmu])-1;
    d = D(tmp_,:);
    d(logical(tmp2_)) = randsample( Nd,sum(tmp2_),'true' );
    D2(ii,:) = d;
end
D = D2;