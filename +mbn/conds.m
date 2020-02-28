function [Cpmnew,B] = conds( Cpm,E,e,B )
%{
Condition a set of cpms (Cpm) by given evidence (E,e)
Input:
    Cpm: Ncpm x 1 cell array of cpm
    E: 1 x Ne array of evidence variable
    e: 1 x Ne array of evidence state
    B: Nx x 1 cell array of basic state indicator
Output:
    Cpmnew: Ncpm x 1 cell array of cpm after conditioning
    B: Nx x 1 cell array of basic state indicator

Ex:
Cpm{1} = cpmjoint(1,(1:2)',[.7 .3]'); B{1} = sparse(2,2); B{1}(1,1) = 1; B{1}(2,2) = 1;
Cpm{2,1} = cpmjoint(2,(1:2)',[.4 .6]'); B{2,1} = sparse(3,2); B{2}(1,1) = 1; B{2}(2,2) = 1; B{2}(3,:) = 1;
Cpm{3} = cpmcond( [3 1 2],1,[1 1 3;2 1 3;1 2 1;2 2 1;1 2 2],[.4 .6 .3 .7 1]' ); B{3} = sparse(2,2); B{3}(1,1) = 1; B{3}(2,2) = 1;
E = [3 2]; e = [1 2];
[Cpmnew,B] = conds( Cpm,E,e,B );
%}

if length(E)~=length(e)
    error('Lengths of E and e must agree')
end

import mbn.*

Cpmnew = cell(size(Cpm));
for ii = 1:length(Cpm)
    [Cpmnew{ii},B] = cond1(Cpm{ii},E,e,B);
end