function X2Y = Mapping( Cx,Cy,B )
%{
Mapping between two event matrices
Ni: # of r.v.'s of interest

Input:
    Cx: Ncx x Ni array of event matrix
    Cy: Ncy x Ni array of event matrix
    B: Ni x 1 cell array of basic event indicator
Output:
    X2Y: Ncx x Ncy logical array of compatibility indicator
%}

import mbn.*
Ncx = size(Cx,1); Ncy = size(Cy,1);

[Cy2,~,ib] = unique(Cy,'rows');
Ncy2 = size(Cy2,1);
X2Y2 = false(Ncx,Ncy2);
for ii = 1:Ncy2
    ci = Cy2(ii,:);
    X2Y2(:,ii) = compats( Cx,ci,B );
end

X2Y = false(Ncx,Ncy);
for ii = 1:Ncy
    X2Y(:,ii) = X2Y2(:,ib(ii));
end

