function [decOpt_prox,decOpt_prox_EVprox_EVcost] = optimize( proxEVsys,EVcost,criticalWeight,Nlay )

Ncomp = length( proxEVsys );
Ndec = cellfun( @length,EVcost );

decOpt_prox = [];
decOpt_prox_EVprox_EVcost = [];

criticalWeight = [0; criticalWeight; criticalWeight(end)+1];
criticalWeight = ( criticalWeight(1:end-1) + criticalWeight(2:end) ) * .5;
for ww = 1:length(criticalWeight)
    weight_l = criticalWeight( ww );
    
    EVsum_sort_w = {};
    EVsum_sortId_w = {};
    for nn = 1:Ncomp
        EVsum_lam_n = EVcost{nn} + weight_l*proxEVsys{nn};
        [EVsum_lam_sort_n,EVsum_lam_sortId_n] = sort( EVsum_lam_n );
        EVsum_sort_w = [EVsum_sort_w; {EVsum_lam_sort_n}];
        EVsum_sortId_w = [EVsum_sortId_w; {EVsum_lam_sortId_n}];
    end

    decOpt_prox_ = zeros( 1,Ncomp );
    EVcost_nextMinRank = ones( 1,Ncomp );
    EVcost_nextMinVal = zeros( 1,Ncomp );
    decOpt_prox_EVprox_ = zeros( 1,Ncomp );
    decOpt_prox_EVcost_ = zeros( 1,Ncomp );

    for nn = 1:Ncomp
        [~,decOpt_sortId] = min( EVsum_sort_w{nn} );
        decOpt_prox_(nn) = EVsum_sortId_w{nn}( decOpt_sortId );
        decOpt_prox_EVprox_(nn) = proxEVsys{nn}( decOpt_prox_(nn) );
        decOpt_prox_EVcost_(nn) = EVcost{nn}( decOpt_prox_(nn) );
        if Nlay > 1
            EVcost_nextMinRank(nn) = EVcost_nextMinRank(nn)+1;
            EVcost_nextMinVal(nn) = EVsum_sort_w{nn}( EVcost_nextMinRank(nn) );
        end
    end
    decOpt_prox = [decOpt_prox; decOpt_prox_];
    decOpt_prox_EVprox_EVcost = [decOpt_prox_EVprox_EVcost; sum(decOpt_prox_EVprox_) sum(decOpt_prox_EVcost_)];

    for oo = 2:Nlay
        [~,nextOptID] = min( EVcost_nextMinVal );
        decOpt_prox_( nextOptID ) = EVsum_sortId_w{nextOptID}( EVcost_nextMinRank(nn) );
        decOpt_prox_EVprox_( nextOptID ) = proxEVsys{nn}( decOpt_prox_( nextOptID ) );
        decOpt_prox_EVcost_( nextOptID ) = EVcost{nn}( decOpt_prox_( nextOptID ) );

        if EVcost_nextMinRank( nextOptID ) == Ndec(nn)
            EVcost_nextMinRank( nextOptID ) = -1;
        else
            EVcost_nextMinRank( nextOptID ) = EVcost_nextMinRank( nextOptID )+1;
            EVcost_nextMinVal( nextOptID ) = EVsum_sort_w{nn}( EVcost_nextMinRank(nextOptID) );
        end
        EVcost_nextMinVal( EVcost_nextMinRank<0 ) = max( EVcost_nextMinVal );

        decOpt_prox = [decOpt_prox; decOpt_prox_];
        decOpt_prox_EVprox_EVcost = [decOpt_prox_EVprox_EVcost; sum(decOpt_prox_EVprox_) sum(decOpt_prox_EVcost_)];
    end
end

[~,decOpt_prox_sortId] = unique( decOpt_prox_EVprox_EVcost,'rows' );
decOpt_prox = decOpt_prox(decOpt_prox_sortId,:);
decOpt_prox_EVprox_EVcost = decOpt_prox_EVprox_EVcost(decOpt_prox_sortId,:);