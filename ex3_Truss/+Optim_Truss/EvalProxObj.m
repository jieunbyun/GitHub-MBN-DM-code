function [EVs_prox, EC, EVs_db] = EvalProxObj(Px,Pl,Cd,db)
%{
Compute the measure required for proposed strategy of optimization

INPUT:
    Px: probability matrix for comp's (members)
    Pl: probability vector for Load
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
Nload = length(Pl); % # of loads being considered
Nfail = size(Px,1)/Nd/Nload; % # of failure modes
EVs_prox = zeros(Nd,Nx);

Nfac = Nfail*Nload; % # of factors
EVs_db = zeros(Nfac,1);

for ff = 1:Nfac
   fl = ceil(ff/Nload);
   ll = ff - Nload*(fl-1);
   
   pb_ = Px( sub2ind( size(Px), Nd*Nload*(fl-1)+Nd*(ll-1)+db, 1:Nx ) );
   evs_ = sum( log(pb_) ) + log(Pl(ll));
   evs_ = exp(evs_);
   if evs_
       delp_ = Px( sub2ind( size(Px), repmat(Nd*Nload*(fl-1)+Nd*(ll-1)+(1:Nd)',1,Nx), repmat( 1:Nx,Nd,1 ) ) );
       delp_ = delp_ ./ repmat( pb_,Nd,1 );   
       EVs_prox = EVs_prox+evs_*delp_;
   end   
   EVs_db(ff) = evs_;
end

EC = Cd - repmat( Cd( sub2ind(size(Cd),db,1:Nx) ),Nd,1 );