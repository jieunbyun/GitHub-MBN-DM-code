function db = RandDec( Nx,Nd )
%{
get decision rule by random generation

INPUT:
    Nx: # of comp's
    Nd: # of decision alternatives
OUTPUT:
    db: 1xNx vector with selected decision alternatives
%}

db = randsample( Nd,Nx,'true' );
db = db(:)';