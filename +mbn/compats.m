function [compcheck,Cred,Bx] = compats( Cx,cx,Bx )
%{
Compatibility check b/w
a set of assingments (Cx) and an assignment (cx)
over multiple r.v.'s

Input:
    Cx: Nr x Nx array of assignments
    cx: 1 x Nx array of assignment
    Bx: Nx x 1 cell array of basic state indicator
Output:
    compcheck: Nr x 1 array of logical index
    Cred: Nr' x 1 array of conditioned assignments
    Bx: Nx x 1 cell array of basic state indicator
%}

import mbn.*

Nr = size(Cx,1); Nxe = length(cx);
Cid = 1:Nr; Cred = zeros(Nr,Nxe);

for ii = 1:Nxe
   [cm,cr,Bi] = compat1(Cx(Cid,ii),cx(ii),Bx{ii});
   Bx{ii} = Bi;
   Cred( cm,ii ) = cr;
   Cred = Cred( cm,: );
   Cid = Cid( cm );
end

compcheck = zeros( Nr,1 ); compcheck(Cid) = 1; compcheck = logical( compcheck );