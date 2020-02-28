%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to obtain probability(Pf) for given magnitude(M).             % 
% Pf is the cell with the same entities as the number of briges.         %
% A cell contains 2 by 4 matrix. In the entity (i,j), i denotes the      % 
% bridge state(1:pristine, 2:damaged), j denotes the damage state        %
% (1:insignificant, 2:moderate, 3:heavy, 4:complete)                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Pf=M2Pi_Rev1_Ac(R_SEIS,M,Ac)

import SF_fragility.*

format long e
% Parameters defined as in Campbell (1997)
% M = from 6.0 to 8.5

F = 0;
S_SR = 0;
S_HR = 0;
D = 1;
FvD = 0;

% Ratio = Sa/PGA from Chopra page 220 
ratio=0.55; 
% ratio=[1.1, 0.75];

% 1 = single bent bridge, 2 =  double bent bridge
% ID=[1, 2, 1, 1, 1, 2, 2, 1, 1, 1];

Pf=cell(length(R_SEIS),1);
Sa=zeros(length(R_SEIS),1);
PGV=zeros(length(R_SEIS),1);
for i=1:length(R_SEIS);
    % Attenuation Relationship by Campbell (1997) 
    mean_ln_PGA(i) = -3.512+0.904*M-1.328*log( sqrt( R_SEIS(i)^2+(0.149*exp(0.647*M))^2 ) )+...
                     (1.125-0.112*log(R_SEIS(i))-0.0957)*F+...
                     (0.440-0.171*log(R_SEIS(i)))*S_SR+...
                     (0.405-0.222*log(R_SEIS(i)))*S_HR; % units of g
    % Sa(i) = exp(mean_ln_PGA(i))*ratio(ID(i)); % units of g
    Sa(i) = exp(mean_ln_PGA(i))*ratio; % units of g
    
    if(D>=1) % units of km 
        FvD = 0;
    else
        FvD = -0.30*(1-S_HR)*(1-D)-0.15*(1-D)*S_SR;
    end
    
    mean_ln_PGV(i) = mean_ln_PGA(i)+0.26+0.29*M...
                     -1.44*log(R_SEIS(i)+0.0203*exp(0.958*M))...
                     +1.89*log(R_SEIS(i)+0.361*exp(0.576*M))...
                     +(0.0001-0.000565*M)*R_SEIS(i)-0.12*F...
                     -0.15*S_SR-0.30*S_SR...
                     +0.75*tanh(0.51*D)*(1-S_HR)+FvD; % units of cm/s
    PGV(i) = exp(mean_ln_PGV(i)); % units of cm/s
    
    Pex_ins = Ind_Fragility_Ac(Sa(i),PGV(i),0.02,Ac);
    Pex_mod = Ind_Fragility_Ac(Sa(i),PGV(i),0.04,Ac);
    Pex_com(1) = Pex_ins(2); Pex_com(2)=0; 
    
    Pf{i}(1,1) = 1-Pex_ins(1); 
    Pf{i}(1,2) = Pex_ins(1)-Pex_mod(1); 
    Pf{i}(1,3) = Pex_mod(1)-Pex_com(1); 
    Pf{i}(1,4) = Pex_com(1);            
end
end