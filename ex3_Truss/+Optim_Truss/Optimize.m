function Dopt = Optimize( EVs_prox,EC,Nlay,D )
%{
Optimization given Evs_prox,EC

INPUT:
    EVs_prox: proxy measure for E[Vs|d] for each decision rule d
              Nd x Nx matrix
    EC: cost E[C|d] (exact)
        Nd x Nx matrix
    Nlay: # of layers for optimal evaluation (To compensate the loss of solutions by weighted sum)
    D: optimal solutions that have been obtained so far
       Nsol x Nx matrix

OUPUT:
    Dopt: optimal decision rules according to proxy measures
          Nsol x Nx matrix
%}

Nx = size(D,2);

lamb = [];
for ii = 1:Nx
    ec_ = EC(:,ii); ev_ = EVs_prox(:,ii);
    [ec_,tmp_] = sort(ec_); ev_ = ev_(tmp_);
        
    ec2_ = []; ev2_ = [];
    while ~isempty( ec_ )
        c_ = ec_(1); v_ = ev_(1);
        if ~all( ec_(2:end)>c_ | ev_(2:end)>v_ )
            ec_(1)=[]; ev_(1)=[];
        else
            ec2_=[ec2_;c_]; ev2_=[ev2_;v_];
            tmp_ = find( ~(ev_(2:end)<v_ | ec_(2:end)<c_) ) + 1;
            ec_( [1 tmp_(:)'] ) = []; ev_( [1 tmp_(:)'] ) = [];
        end
    end
        
    lamb2_ = [];
    while length( ec2_ ) > 1
        lam_ = (ec2_(1) - ec2_(2:end)) ./ (ev2_(1) - ev2_(2:end));
        [lam_,tmp_] = max( lam_ );
        lamb2_ = [lamb2_ lam_]; 
        ec2_(1:tmp_) = []; ev2_(1:tmp_) = [];
    end
    
    if ~isempty( lamb2_)
%         lamb2_ = -1./lamb2_;
        lamb2_ = -lamb2_;
        lamb2_ = ([0 lamb2_] + [lamb2_ lamb2_(end)+1])*.5;

        lamb = [lamb lamb2_];
    end
end

lamb = unique(lamb);
lamb = sort(lamb);
Nlam = length(lamb);

Dopt = [];
for ii = 1:Nlam
   E_ = EC + lamb(ii)*EVs_prox;
   [~,dopt_] = min(E_);
   D2 = dopt_;
   E_( sub2ind(size(E_),dopt_,1:Nx) ) = max( E_(:) );
   for jj = 1:Nlay        
       dopt2_ = dopt_; % assume that Nlay < Nx
       [r,c] = FindMinIdx(E_);
       dopt2_(c) = r;       
       D2 = [D2; dopt2_];       
       E_( sub2ind(size(E_),r,c) ) = max( E_(:) );
   end
   tmp_ = ~ismember(D2,[Dopt; D],'rows');
   Dopt = [Dopt;D2(tmp_,:)];
end

function [r,c] = FindMinIdx(E)
%{
find indices of row and column for matrix E's min value

Input: E: matrix
Output: r,c: row and column indices of matrix E's min value
%}
[~,I] = min(E(:));
[r,c] = ind2sub(size(E),I);