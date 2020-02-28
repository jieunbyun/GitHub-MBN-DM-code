function [Ncomp, Cpm, var, B, val, EVcost] = MBNquant_RBD


import mbn.*
Ncomp = 8;

var_ = 0;
% Dn
Ndec = 3;
eventMat_D = (1:Ndec)'; probVec_D = ones(Ndec,1); val_D = {'comp1' 'comp2' 'comp3'}'; B_D = eye(Ndec);
for nn = 1:Ncomp
   var_ = var_+1;
   var.D(nn) = var_;
   Cpm{var_,1} = cpmjoint( var_,eventMat_D,probVec_D );
   B{var_,1} = B_D;
   val{var_,1} = val_D;
end

% Xn
compType1 = [1 2 4 7];
eventMat_X = [repmat((1:2)',Ndec,1) repelem(eventMat_D,2,1)];
probVec_X1 = [.99 .01 .997 .003 .999 .001]';
probVec_X2 = [.99 .01 .995 0.005 .999 .001]';
val_X = {'Survive' 'Fail'}'; B_X = [eye(2); 1 1];
for nn = 1:Ncomp
   var_ = var_+1;
   var.X(nn) = var_;
   if ismember( nn,compType1 )
       probVec_X_n = probVec_X1;
   else
       probVec_X_n = probVec_X2;
   end
       
   Cpm{var_,1} = cpmcond( [var_ var.D(nn)],1,eventMat_X,probVec_X_n );
   B{var_,1} = B_X;
   val{var_,1} = val_X;
end

% S: X_(N+1)
eventMat_S = 3*ones( 5,1+Ncomp ); 
eventMat_S( 1,[1 9] ) = [2 2]; eventMat_S( 2,[1 8 9] ) = [2 2 1]; eventMat_S( 3,[1 2 3 4 5 8 9] ) = [2 2 2 2 2 1 1];
eventMat_S( 4,[1 2 3 4 5 6 8 9] ) = [2 2 2 2 1 2 1 1]; eventMat_S( 5,[1 2 3 4 5 6 7 8 9] ) = [2 2 2 2 1 1 2 1 1];
probVec_S = ones( size(eventMat_S,1),1 );
val_S = {'Survive' 'Fail'}'; B_S = eye(2);
var_ = var_+1;
var.S = var_;
Cpm{var_,1} = cpmcond( [var_ var.X],1,eventMat_S,probVec_S );
B{var_,1} = B_S;
val{var_,1} = val_S;

% Cost
EVcost1 = [100 150 250]'; EVcost2 = [80 140 200]';
EVcost = {};
for nn = 1:Ncomp
    if ismember(nn,compType1)
        EVcost = [EVcost; {EVcost1}];
    else
        EVcost = [EVcost; {EVcost2}];
    end
end
 