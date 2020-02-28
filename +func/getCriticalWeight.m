function criticalWeight = getCriticalWeight( proxEVsys,EVcost )

Ncomp = length( proxEVsys );

criticalWeight = [];
for nn = 1:Ncomp
   proxEVsys_n =  proxEVsys{nn}; EVcost_n = EVcost{nn};
   [proxEVsys_n,proxEVsys_n_sortId] = sort( proxEVsys_n ); EVcost_n = EVcost_n( proxEVsys_n_sortId );
   
   lam_n = [];
   id_ = 1;
   while id_ < length( proxEVsys_n )
      lam_id_ =  ( EVcost_n( (id_+1):end ) - EVcost_n( id_ ) ) ./ ( proxEVsys_n( (id_+1):end ) - proxEVsys_n( id_ ) );
      [lam_id_min_,lam_id_min_id_] = min( lam_id_ );
      if lam_id_min_ < 0
          lam_n = [lam_n; lam_id_min_];
          id_ = lam_id_min_id_ + id_;
      else
          id_ = length( proxEVsys_n );
      end
   end
   
   criticalWeight = [criticalWeight; -lam_n];   
end
criticalWeight = sort(criticalWeight);

