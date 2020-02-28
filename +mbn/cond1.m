function [cpmnew,B] = cond1(cpm,E,e,B)
%{
Condition a cpm (cpm) by given evidence (E,e)
Input:
    cpm: a cpm
    E: 1 x Ne array of evidence variable
    e: 1 x Ne array of evidence state
    B: Nx x 1 cell array of basic state indicator
Output:
    cpmnew: a cpm after conditioning
    B: Nx x 1 cell array of basic state indicator

Ex:
cpm = cpmcond( [3 1 2],1,[1 1 3;2 1 3;1 2 1;2 2 1;1 2 2],[.4 .6 .3 .7 1]' );
B{1} = sparse(2,2); B{1}(1,1) = 1; B{1}(2,2) = 1;
B{2,1} = sparse(3,2); B{2}(1,1) = 1; B{2}(2,2) = 1; B{2}(3,:) = 1;
B{3} = sparse(2,2); B{3}(1,1) = 1; B{3}(2,2) = 1;
E = [3 2]; e = [1 2];
[cpmnew,B] = cond1( cpm,E,e,B );
%}

if length(E)~=length(e)
    error('Lengths of E and e must agree')
end

import mbn.*

[xe,ic,ie] = intersect(cpm.scope,E);
Nr = size(cpm.C,1); Nxe = length(xe);
Cid = 1:Nr; Cnew2 = zeros(Nr,Nxe);
for ii = 1:Nxe
   xi = xe(ii); 
   [cm,cr,Bi] = compat1(cpm.C(Cid,ic(ii)),e(ie(ii)),B{xi});
   B{xi} = Bi;
   Cnew2( cm,ii ) = cr;
   Cnew2 = Cnew2( cm,: );
   Cid = Cid( cm );
end

Nx = size(cpm.C,2);
Cnew = zeros(length(Cid),Nx);
icn = setdiff( 1:Nx,ic );
Cnew(:,icn) = cpm.C(Cid,icn);
Cnew(:,ic) = Cnew2;

cpmnew = cpm;
cpmnew.C = Cnew;
cpmnew.p = cpm.p(Cid);