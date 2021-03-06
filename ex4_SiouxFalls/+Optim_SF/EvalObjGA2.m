function [EVs, EC] = EvalObjGA2( D,S,S_,Path,Px,Pm,Pl,Cd )
%{
Evaluate solutions (upper bound of P(s0))
INPUT:
    D: decision rules to be evaluated - Nsol x Nx matrix
    S: survival nodes identified by RDA = Nset x 1 cell
    S_: failed nodes identified by RDA - Nset x 1 cell
    Path: path of each link-set
    Px: probability matrix for comp's (failure only)
    Pm: probability vector for magnitude of EQ
    Pl: probability vector for location of EQ
    Cd: cost matrix for comp's
OUTPUT:
    EVs: E[Vs|d] - Nsol x 1 vector 
    EC: E[C|d] - Nsol x 1 vector
%}

Nsol = size(D,1);
Nx = size(D,2);
Nm = length(Pm);
Nl = length(Pl);
Nd = size(Cd,1);
Nset = length(S); % # of failure modes
Nml = Nm*Nl;
Nfac = Nset*Nml;

EC = zeros(Nsol,1);
EVs = zeros(Nsol,1);

for ii = 1:Nsol
   d_ = D(ii,:);
   evs_ = 0;
   for ff = 1:Nfac
      fl = ceil(ff/Nml);
      ff_ = ff - Nml*(fl-1);
      mm = ceil(ff_/Nl);
      ll = ff_ - Nl*(mm-1);
      
      c_ = [S{fl}(:)' Path{fl}(:)' S_{fl}(:)']; % survived comp's failed comp's
      Ns_ = length(S{fl}) + length( Path{fl} );

      p_ = Px( sub2ind( size(Px), Nd*Nl*(mm-1)+Nd*(ll-1)+d_(c_), c_) );
      p_(1:Ns_) = 1-p_(1:Ns_);
      evs_ = evs_+ exp( sum( log(p_)) + log( Pm(mm) ) + log(Pl(ll)));   
   end
   EVs(ii) = 1-evs_;
   EC(ii) = sum( Cd( sub2ind(size(Cd),d_,1:Nx) ) );
end
