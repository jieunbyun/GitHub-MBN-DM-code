%{
Genetic algorithm for Truss bridge example
%}

clear; close all;
import Optim_Truss.*
load Quant_Truss
load Optim_Truss
rng(1)
%%
Npop = 50; % # of populations at each gen.
Ngen = 1e3; % # of gen's to be examined
Rmu = .4; % Ratio of mutated populations
Rcr = .4; % Ratio of pupoluations being cross-over
rmu = .15; % Ratio of mutation (genetic changes)

Nd = size(Cd,1); Nx = size(D,2);
Nmu = floor( Npop*Rmu ); Ncr = floor( Npop*Rcr ); Nran = Npop-Nmu-Ncr;

Nsol = zeros(Ngen,1); 

% 1st gen.
Dga = GenRand( Npop,Nd,Nx );
[EVsga, ECga] = EvalObjGA( Dga,Px,Pl,Cd );
[Dga,EVsga,ECga] = SortNonDominSol( Dga,EVsga,ECga );
Nsol(1) = size(Dga,1);

% later gen.
for ii = 2:Ngen
    D2 = GenRand( Nran,Nd,Nx );
    D2 = [D2;GenMut( Nmu,Nd,Dga,rmu ) ];
    D2 = [D2;GenCro( Ncr,Dga )];
    
    [EVs2, EC2] = EvalObjGA( D2,Px,Pl,Cd );
    [Dga,EVsga,ECga] = SortNonDominSol( [Dga;D2],[EVsga;EVs2],[ECga;EC2] );
    
    Nsol(ii) = size(Dga,1);
    
    if ~rem(ii,1e2)
        eval( ['save GA_Truss_gen' num2str(ii) ' Dga EVsga ECga'] )
        semilogy( EC,EVs,'*' ); hold on; semilogy( ECga,EVsga,'sq' ); grid on; hold off; drawnow;
    end    
    
end
