function D = GenCro( Ncr,D )
%{
Generate Ncr decision rules by cross-over from the ones in D
INPUT: 
    Ncr: # of decision rules to be generated
    D: a set of rules from which the new generation is produced
OUTPUT:
    D: Ncr x Nx matrix with generated rules
%}

Nsol = size(D,1); Nx = size(D,2); D2 = zeros(Ncr,Nx);
for ii = 1:Ncr
   tmp_ =  randsample(Nsol,2,'true');
   tmp2_ = randsample(2,Nx,'true')-1; tmp2_ = logical(tmp2_);
   d = D(tmp_(1),:);
   d(tmp2_) = D(tmp_(2),tmp2_);
   D2(ii,:) = d;
end
D = D2;