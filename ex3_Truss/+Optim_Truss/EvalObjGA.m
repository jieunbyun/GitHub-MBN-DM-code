function [EVs, EC] = EvalObjGA( D,Px,Pl,Cd )
%{
Evaluate solutions
INPUT:
    D: decision rules to be evaluated - Nsol x Nx matrix
    Px: probability matrix for comp's (members)
    Pl: probability vector for Load
    Cd: cost matrix for comp's
OUTPUT:
    EVs: E[Vs|d] - Nsol x 1 vector
    EC: E[C|d] - Nsol x 1 vector
%}

Nsol = size(D,1);
Nx = size(D,2);
Nload = length(Pl);
Nd = size(Cd,1);
Nfail = size(Px,1)/Nd/Nload; % # of failure modes
Nfac = Nfail*Nload;

EC = zeros(Nsol,1);
EVs = zeros(Nsol,1);

for ii = 1:Nsol
   d_ = D(ii,:);
   evs_ = 0;
   for ff = 1:Nfac
      fl = ceil(ff/Nload);
      ll = ff - Nload*(fl-1);

      p_ = Px( sub2ind( size(Px), Nd*Nload*(fl-1)+Nd*(ll-1)+d_, 1:Nx ) );
      evs_ = evs_+ exp( sum( log(p_)) + log(Pl(ll)));    
   end
   EVs(ii) = evs_;
   EC(ii) = sum( Cd( sub2ind(size(Cd),d_,1:Nx) ) );
end
