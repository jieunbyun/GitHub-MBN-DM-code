function Prob = GetProb(F,Sig,Area,LoadRat)
%{
Get the probabilities of probabilities of members' survival

INPUT:
    F: demand on each members 1xNmem
    Sig: structure of info. on yield stress (mu and cov)
         consistent for all members
    Area: cross-sectional areas to be examined 1xNarea
          consistent for all members
    LoadRat: Scale for loads to be examined w.r.t. the load for which F has been evaluated
             (LoadRat(1) = 1) 1xNload
    
OUTPUT:
    Prob: (Narea x Nload) x Nmem matrix for
            (Load Dn) = [1 1;1 2;...;1 Narea;2 1;2 2;...;2 Narea;...;Nload 1;Nload 2;...;Nload Narea]
%}
Narea = length(Area);
Nload = length(LoadRat);
Nmem = length(F);
F = abs(F);
F = F(:)'; Area = Area(:); LoadRat = LoadRat(:);

F = repmat(F,Nload,1) .* repmat( LoadRat,1,Nmem );
F = repelem(F,Narea,1);
Area = repmat(Area,Nload,Nmem);
FA = F./Area;
Prob = 1-normcdf( FA,Sig.mu,Sig.cov*Sig.mu );