function [EVs, EC, EVs_prox, AG] = EvalObj( D,S,S_,Px,Pm,Pl,Cd,EVs_prox,EVs_db,db )
%{
Evaluate exact values of E[Vs|d], E[C|d], and measure A/G for selecting
next basis decision rule

INPUT:
    D: decision rules to be evaluated - Nsol x Nx matrix
    S: survival nodes identified by RDA = Ncut x 1 cell
    S_: failed nodes identified by RDA - Ncut x 1 cell
    Px: probability matrix for comp's (failure only)
    Pm: probability vector for magnitude of EQ
    Pl: probability vector for location of EQ
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
Nm = length(Pm);
Nl = length(Pl);
Nd = size(Cd,1);
Ncut = length(S); % # of failure modes
Nml = Nm*Nl;
Nfac = Ncut*Nml;
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
          fl = ceil(ff/Nml);
          ff_ = ff - Nml*(fl-1);
          mm = ceil(ff_/Nl);
          ll = ff_ - Nl*(mm-1);
          
          c_ = [S{fl}(:)' S_{fl}(:)']; % failed comp's / survived comp's
          Ns_ = length(S{fl});

          p_ = Px( sub2ind( size(Px), Nd*Nl*(mm-1)+Nd*(ll-1)+d_(c_), c_ ) );
          p_(1:Ns_) = 1-p_(1:Ns_);
          evs_ = evs_+ exp( sum( log(p_)) + log( Pm(mm) ) + log(Pl(ll)) );

          delpb_ = Px( sub2ind( size(Px), Nd*Nl*(mm-1)+Nd*(ll-1)+db(c_), c_ ) );
          delpb_(:,1:Ns_) = 1-delpb_(:,1:Ns_);
          delp_ = (p_ ./ delpb_-1);
          geo_ = geo_ + EVs_db(ff)*exp( sum( log(delp_+1) ) );
          ari_ = ari_ + EVs_db(ff)*sum(delp_+1);
    

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
          fl = ceil(ff/Nml);
          ff_ = ff - Nml*(fl-1);
          mm = ceil(ff_/Nl);
          ll = ff_ - Nl*(mm-1);
          
          c_ = [S{fl}(:)' S_{fl}(:)']; % failed comp's / survived comp's
          Ns_ = length(S{fl});

          p_ = Px( sub2ind( size(Px), Nd*Nl*(mm-1)+Nd*(ll-1)+d_(c_), c_) );
          p_(1:Ns_) = 1-p_(1:Ns_);
          evs_ = evs_+ exp( sum( log(p_)) + log( Pm(mm) ) + log(Pl(ll)));    
       end
       EVs(ii) = evs_;
       EC(ii) = sum( Cd( sub2ind(size(Cd),d_,1:Nx) ) );
       EVs_prox(ii) = sum( EVs_prox_( sub2ind(size(EVs_prox_),d_(c_),c_) ) );
    end
end