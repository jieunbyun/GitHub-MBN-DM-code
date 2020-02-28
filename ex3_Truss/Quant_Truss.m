%{
2018-11-17 Quantify Px and Cd in Truss example

Px: (Nfail x Narea x Nload) x Nmem matrix for
    (Load Dn) = [1 1;1 2;...;1 Narea;2 1;2 2;...;2 Narea;...;Nload 1;Nload 2;...;Nload Narea]
Cd: Narea x Nmem matrix
    Cost matrix for decision alternatives
Fmode: Nfail x 1 vector for Failure modes
       The order of elements indicate the order of failures of members
       (Sequential failures are of interest)
%}

clear; close all;
import SA_Truss.*
import Optim_Truss.*
%% Structural info.
% Load
Load.mu = 70; % kN
Load.cov = .1; 
Load.Ngrid = 9; % # of grids (odd)
Load.gridw = 5; % width of grid
Load.val = Load.mu-Load.gridw*(Load.Ngrid/2-1)+Load.gridw*(0:(Load.Ngrid-2));
Load.p = normcdf( Load.val,Load.mu,Load.mu*Load.cov );
Load.p = diff( [0 Load.p 1] ); Load.p = Load.p';
Load.val = Load.mu - floor(Load.Ngrid/2)*Load.gridw + Load.gridw*(0:Load.Ngrid-1);
Load.ratio = Load.val ./ Load.val(1);

% Yield stress
Sig.mu = 276e3; % kN/m^2
Sig.cov = .05;

% Cross section
Area.Ngrid = 7;
Area.val = (12+0.5*((1:Area.Ngrid)'-1))*1e-4; % (m^2)

% Structure
D.Con = [1 2;2 3;3 4;4 5;5 6;6 7;1 8;8 9;9 10;10 11;11 12;7 12;2 8;3 9;4 10;5 11;6 12;...
    3 8;2 9;4 9;3 10;4 11;5 10;5 12;6 11]'; % member - node
Nmem = size(D.Con,2); % # of members

D.Coord = [0 0;2 1.6;4 2.2;6 2.6;8 2.2;10 1.6;12 0;2 0;4 0;6 0;8 0;10 0];
Nnod = size(D.Coord,1); % # of nodes
D.Coord = [D.Coord zeros( Nnod,1 )]';

D.Re = zeros( Nnod,3 ); D.Re(1,:) = 1; D.Re(7,2) = 1; D.Re(:,3) = 1; D.Re = D.Re';
D.Load = zeros( Nnod,3 ); D.Load( 8:12,2 ) = -Load.val(1); D.Load = D.Load'; % kN (smallest value)

D.E = ones( 1,Nmem )*1e4; D.A = ones( 1,Nmem ); % dummy (only Force of interest)

%% Qunatify P
F1 = ST(D);
P1 = GetProb(F1,Sig,Area.val,Load.ratio);

Px = [];
Fmode = {}; % failure modes
for ii = 1:Nmem
   D2 = D;
   D2.E(ii) = 0;
   F2 = ST(D2);
   if sum( isnan(F2) )
       Fmode = [Fmode; ii];
       px_ = P1;
       px_(:,ii) = 1-P1(:,ii);
       Px = [Px; px_];
   else
       P2 = GetProb(F2,Sig,Area.val,Load.ratio);
       for jj = setdiff( 1:Nmem,ii )
           D3 = D2;
           D3.E(jj) = 0;
           F3 = ST(D3);
           p2_ = P2(:,jj)-P1(:,jj);
           if ~sum(p2_<0)
               if sum(isnan(F3))
                   px_ = P2;
                   px_(:,ii) = 1-P1(:,ii);
                   px_(:,jj) = p2_;
                   Fmode = [Fmode; [ii jj]];
                   Px = [Px;px_];           
               else
                   P3 = GetProb(F3,Sig,Area.val,Load.ratio);
                   for kk = setdiff( 1:Nmem,[ii jj] )
                       p3_ = P3(:,kk)-P2(:,kk);
                       if ~sum(p3_<0)
                           Fmode = [Fmode; [ii jj kk]];
                           px_ = P3;
                           px_(:,ii) = 1-P1(:,ii);
                           px_(:,jj) = p2_;
                           px_(:,kk) = p3_;
                           Px = [Px;px_];
                       end
                   end
               end
           end
       end
   end
end

%% Quantify C
Carea = [10 30 45 55 60 62.5 65]'; % Cost/unit length (m) for each area
Cd = [];
for ii=1:Nmem
   H=D.Con(:,ii);C=D.Coord(:,H(2))-D.Coord(:,H(1));Le=norm(C);   
   Cd = [Cd Carea.*Le/1e3];
end

Pl = Load.p;
save Quant_Truss Px Fmode Cd Pl