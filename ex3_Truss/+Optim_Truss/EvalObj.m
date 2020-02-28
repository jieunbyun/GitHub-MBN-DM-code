function [EVs, EC, EVs_prox, AG] = EvalObj( D,Px,Pl,Cd,EVs_prox,EVs_db,db )
%{
Evaluate exact values of E[Vs|d], E[C|d], and measure A/G for selecting
next basis decision rule

INPUT:
    D: decision rules to be evaluated - Nsol x Nx matrix
    Px: probability matrix for comp's (members)
    Pl: probability vector for Load
    Cd: cost matrix for comp's
    EVs_prox: Nd x Nx matrix for proxy measure of each Val(Dn) for each Dn
    EVs_db: E[Vs(c)|db] for each factor (exact) - Nfac x 1 vector
    db: basis decision rule - 1 x Nx vector
    
OUTPUT:
    EVs: E[Vs|d] - Nsol x 1 vector
    EC: E[C|d] - Nsol x 1 vector
    EVs_prox: E_tilde[Vs|d] - Nsol x 1 vector
    AG: measure A/G for selecting next basis decision rule - Nsol x 1 vector
%}

Nsol = size(D,1);
Nx = size(D,2);
Nload = length(Pl);
Nd = size(Cd,1);
Nfail = size(Px,1)/Nd/Nload; % # of failure modes
Nfac = Nfail*Nload;
EVs_prox_ = EVs_prox;

EC = zeros(Nsol,1);
EVs = zeros(Nsol,1);
EVs_prox = zeros(Nsol,1);

if nargout>3
    G = zeros(Nsol,1);
    A = zeros(Nsol,1);

    for ii = 1:Nsol
       d_ = D(ii,:);
       geo_ = 0; ari_ = 0;
       evs_ = 0;
       for ff = 1:Nfac
          fl = ceil(ff/Nload);
          ll = ff - Nload*(fl-1);

          p_ = Px( sub2ind( size(Px), Nd*Nload*(fl-1)+Nd*(ll-1)+d_, 1:Nx ) );
          evs_ = evs_+ exp( sum( log(p_)) + log(Pl(ll)));

          if EVs_db(ff)
              delpb_ = Px( sub2ind( size(Px), Nd*Nload*(fl-1)+Nd*(ll-1)+db, 1:Nx ) );
              delp_ = (p_ ./ delpb_-1);
              geo_ = geo_ + EVs_db(ff)*exp( sum( log(delp_+1) ) );
              ari_ = ari_ + EVs_db(ff)*sum(delp_+1);
          end      

       end
       EVs(ii) = evs_;
       G(ii) = geo_; A(ii) = ari_;
       EC(ii) = sum( Cd( sub2ind(size(Cd),d_,1:Nx) ) );
       EVs_prox(ii) = sum( EVs_prox_( sub2ind(size(EVs_prox_),d_,1:Nx) ) );
    end
    AG = A./G;
else
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
       EVs_prox(ii) = sum( EVs_prox_( sub2ind(size(EVs_prox_),d_,1:Nx) ) );
    end
end