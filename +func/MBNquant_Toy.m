function [Ncomp, Cpm, var, B, val, EVcost] = MBNquant_Toy

import mbn.*
Ncomp = 3;

var_ = 0;
% Dn
eventMat_D = (1:2)'; probVec_D = [1 1]'; val_D = {'Nothing' 'Retrofit'}'; B_D = eye(2);
for nn = 1:Ncomp
   var_ = var_+1;
   var.D(nn) = var_;
   Cpm{var_,1} = cpmjoint( var_,eventMat_D,probVec_D );
   B{var_,1} = B_D;
   val{var_,1} = val_D;
end

% H
var_ = var_+1;
eventMat_H = (1:2)'; probVec_H = [.2 .8]'; val_H = {'Occur' 'None'}'; B_H = eye(2);
var.H = var_;
Cpm{var_,1} = cpmjoint( var_,eventMat_H,probVec_H );
B{var_,1} = B_H;
val{var_,1} = val_H;

% Xn
eventMat_X = [repmat(eventMat_D,2,1) repelem(eventMat_H,2,1)];
eventMat_X = [repmat((1:2)',2*2,1) repelem(eventMat_X,2,1)];
probVec_X = [.6 .4 .7 .3 .8 .2 .9 .1]';
val_X = {'Survive' 'Fail'}'; B_X = [eye(2); 1 1];
for nn = 1:Ncomp
   var_ = var_+1;
   var.X(nn) = var_;
   Cpm{var_,1} = cpmcond( [var_ var.D(nn) var.H],1,eventMat_X,probVec_X );
   B{var_,1} = B_X;
   val{var_,1} = val_X;
end

% S: X_(N+1)
eventMat_S = [2 2 2 3; 2 2 1 2];
probVec_S = [1 1]';
val_S = {'Survive' 'Fail'}'; B_S = eye(2);
var_ = var_+1;
var.S = var_;
Cpm{var_,1} = cpmcond( [var_ var.X],1,eventMat_S,probVec_S );
B{var_,1} = B_S;
val{var_,1} = val_S;

% Cost
EVcost = {[0 100]'};
EVcost = [EVcost; {[0 60]'}];
EVcost = [EVcost; {[0 50]'}];
