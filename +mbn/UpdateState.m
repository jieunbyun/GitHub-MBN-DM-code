function [B,id] = UpdateState( b,B )
    if ismember( b,B,'rows' )
        id = find( ismember( B,b,'rows' ) );
    else
        id = size(B,1)+1;
        B = [B;b];
    end
end