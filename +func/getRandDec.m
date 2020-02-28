function randDec = getRandDec( Ndec )

Ncomp = length(Ndec);
randDec = zeros( 1,Ncomp );
for nn = 1:length(Ndec)
    randDec(nn) = randsample( Ndec(nn),1 );
end