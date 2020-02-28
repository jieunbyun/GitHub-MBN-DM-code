function D = GenRand( Npop,Nd,Nx )
%{
Generate Npop decision rules randomly
INPUT: 
    Npop: # of decision rules to be generated
    Nd: # of decision alternatives
    Nx: # of decision variables
OUTPUT:
    D: Npop x Nx matrix with generated rules
%}

D = randsample(Nd,Npop*Nx,'true');
D = reshape(D,Npop,Nx);