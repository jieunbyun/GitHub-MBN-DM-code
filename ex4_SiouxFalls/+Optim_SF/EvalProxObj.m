function [EVs_prox, EC, EVs_db] = EvalProxObj(S,S_,Px,Pm,Pl,Cd,db)
%{
Compute the measure required for proposed strategy of optimization

INPUT:
    S: survival nodes of cuts identified by RDA = Ncut x 1 cell
    S_: failed nodes of cuts identified by RDA - Ncut x 1 cell
    Px: probability matrix for comp's (failure only)
    Pm: probability vector for magnitude of EQ
    Pl: probability vector for location of EQ
    Cd: cost matrix for comp's
    db: basis decision rule
OUTPUT:
    EVs_prox: proxy measure for E[Vs|d] for each decision rule d
              Nd x Nx matrix
    EC: cost E[C|d] (exact)
        Nd x Nx matrix
    EVs_db: E[Vs(c)|db] for each factor (exact)
            Nfac x 1 vector
    
%}
Nd = size(Cd,1); % # of decisions
Nx = size(Cd,2); % # of comp's 
Nm = length(Pm);
Nl = length(Pl);
Ncut = length(S); % # of cut-sets

Nml = Nm*Nl;
Nfac = Ncut*Nml;
EVs_prox = zeros(Nd,Nx);

EVs_db = zeros(Nfac,1);

for ff = 1:Nfac
   fl = ceil(ff/Nml);
   ff_ = ff - Nml*(fl-1);
   mm = ceil(ff_/Nl);
   ll = ff_ - Nl*(mm-1);
   
   c_ = [S{fl}(:)' S_{fl}(:)']; % failed comp's / survived comp's
   Ns_ = length(S{fl}); Nc_ = length(c_);
   
   pb_ = Px( sub2ind( size(Px), Nd*Nl*(mm-1)+Nd*(ll-1)+db(c_), c_ ) );
   pb_(1:Ns_) = 1-pb_(1:Ns_);
   delp_ = Px( sub2ind( size(Px), repmat(Nd*Nl*(mm-1)+Nd*(ll-1)+(1:Nd)',1,Nc_), repmat( c_,Nd,1 ) ) );
   delp_(:,1:Ns_) = 1-delp_(:,1:Ns_);
   delp_ = delp_ ./ repmat( pb_,Nd,1 );  
      
   evs_ = log( Pm(mm) ) + log( Pl(ll) ) +sum( log(pb_) );   
   evs_ = exp(evs_);
   
   EVs_prox( :,c_ ) = EVs_prox( :,c_ ) + evs_*delp_;
   EVs_db(ff) = evs_;
end

EC = Cd - repmat( Cd( sub2ind(size(Cd),db,1:Nx) ),Nd,1 );