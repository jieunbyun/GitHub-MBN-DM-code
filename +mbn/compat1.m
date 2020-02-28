function [compcheck,Cred,Bx] = compat1( Cx,cx,Bx )
%{
Compatibility check b/w
a set of assingments (Cx) and an assignment (Cx)
over a single r.v.

Input:
    Cx: Nr x 1 array of assignments
    cx: scalar of assignment
    Bx: Ns x Nbs array of basic state indicator
Output:
    compcheck: Nr x 1 array of logical index
    Cred: Nr' x 1 array of conditioned assignments
    Bx: Ns' x Nbs array of basic state indicator

Ex:
Cx = (1:7)'; cx = 5;
Bx = sparse(7,4); Bx(1,1) = 1; Bx(2,2) = 1; Bx(3,3) = 1; Bx(4,4) = 1; Bx(5,1:3) = 1; Bx(6,2:4) = 1; Bx(7,:) = 1;
[compcheck,Cred] = compat1( Cx,cx,Bx );
%}

import mbn.*
Nb = size(Bx,2); Ns = size(Bx,1); Nc = length(Cx);
compcheck = zeros(Nc,1);
bc = Bx(cx,:); Cred = [];
for ii = 1:Nc
   ei = bc.*Bx(Cx(ii),:);
   if sum(ei)
       compcheck(ii) = 1;
       ci = zeros(1,Nb);
       ci(ei>0) = 1;
       ei2 = find( ismember( Bx,ci,'rows' ) );
       if isempty(ei2)
           Ns = Ns+1;
           ei2 = Ns;
           Bx = [Bx;sparse(ci)];
       end
       Cred = [Cred;ei2];
   end   
end
compcheck = logical(compcheck);