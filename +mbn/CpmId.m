function cid = CpmId( xid,Cpm )
%{
find cpm index having xid
%}

import mbn.*
cid = [];
for ii = 1:length(Cpm)
    cpm_ = Cpm{ii};
    if isa(cpm_,'cpmjoint')
        if ismember(xid,cpm_.scope)
            cid = [cid ii];
        end
    else
        if ismember(xid,cpm_.scope(1:cpm_.nc))
            cid = [cid ii];
        end
    end
end