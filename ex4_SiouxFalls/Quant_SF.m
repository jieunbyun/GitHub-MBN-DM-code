%{
2018-11-17 Quantify Px, Pm, Pl, Cd in Sioux Falls Benchmark Network

Px: (Nfail x Narea x Nload) x Nmem matrix for
    (Load Dn) = [1 1;1 2;...;1 Narea;2 1;2 2;...;2 Narea;...;Nload 1;Nload 2;...;Nload Narea]
Cd: Narea x Nmem matrix
    Cost matrix for decision alternatives
Fmode: Nfail x 1 vector for Failure modes
       The order of elements indicate the order of failures of members
       (Sequential failures are of interest)
%}

clear; close all;
import SF_fragility.*
import Optim_SF.*
%% Pm, Pl, Cd
Nx = 76;
% Magnitude of EQ
Mag.val = linspace(6,8.5,6);
Mag.p = 1 / ( 1-exp( -0.76*(8.5-6) ) ) * (1-exp(-0.76*(Mag.val-6)));
Mag.p = diff(Mag.p(:));
Mag.val = ( Mag.val(1:end-1) + Mag.val(2:end) )' * .5;
Nm = length(Mag.val);

% Location of EQ
Loc.val = 0.5*[(-4:5)' (5:-1:-4)']; % EQ location (x,y)
Nl = length(Loc.val);
Loc.p = ones( Nl, 1 ) / Nl; % Uniform dist

% Target Capacity (PGA)
lnPga = log( 0.8+(0:0.1:0.3)' ); % log(Pga) - unit:g
Nd = length(lnPga);

% Cost
Cost = [50 80 120 170]'; % Cost for each area
Cd = repmat(Cost,1,Nx);

%% Qunatify Px
load SiouxFalls_RDA G
SF.coord = load('data/SiouxFalls_node.txt');
SF.coord = SF.coord(:,2:3)*0.000025; % unit :in -> km 
SF.arc = load('data/SiouxFalls_net1.txt');
SF.coord = ( SF.coord( SF.arc(:,1),: ) + SF.coord( SF.arc(:,2),: ) ) * .5; % Coord of arcs regarded as mid-point

Px = zeros( Nm*Nl*Nd,Nx ); % Only fail prob.
for dd = 1:Nd
    d_ = lnPga(dd);
    for ll = 1:Nl
        l_ = Loc.val(ll,:);
        r_seis = sqrt( ( SF.coord(:,1)-l_(1) ).^2 + ( SF.coord(:,2)-l_(2) ).^2 );
        for mm = 1:Nm
            m_ = Mag.val(mm);
            %Attenuation law: Campbell (1997)
            lnD_ = -3.512 + 0.904 * m_ -1.328 * log( sqrt(r_seis.^2 + (0.149*exp(0.647*m_))^2) ); % median demand (pga,g)
            sig_ = 0.173 - 0.140*lnD_; % standard error of ln(Ah)
            sig_( exp(lnD_)<.068 ) = 0.55; sig_( exp(lnD_)>0.21 ) = 0.39; 
            Pf_ = normcdf( (lnD_-d_)./sig_ );
            tmp_ = Nl*Nd*(mm-1)+Nd*(ll-1)+dd;
            Px( sub2ind( size(Px),tmp_*ones(1,Nx),1:Nx ) ) = Pf_;
        end
    end
end

%% save data
Pm = Mag.p; Pl = Loc.p;

tmp = cellfun( @isempty,G.PATHNODE ); % cut-sets
S = G.S(tmp); S_ = G.S_(tmp);
S2 = G.S(~tmp); S2_ = G.S_(~tmp); Path2 = G.PATHLINK(~tmp);  % link-sets
save Quant_SF Px Pm Pl Cd S S_
save Quant_SF2 Px Pm Pl Cd S2 S2_ Path2