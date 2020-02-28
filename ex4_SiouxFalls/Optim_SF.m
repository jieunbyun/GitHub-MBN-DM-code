%{
2018-11-20 Optimize ID for SF example (using cut-sets)

D: optimal decision rules - Nopt x Nx matrix
EVs: system failure probability - Nopt x 1 vector
EC: cost - Nopt x 1 vector
Db: utilized basis decision rules - Ndb x Nx matrix
%}

clear; close all;
import Optim_SF.*
load Quant_SF

rng(1)
%%
Ndb = 50; % # of basis decision rules
Nx = size(Cd,2); % # of comp's
Nd = size(Cd,1); % # of decision alternatives
Nlay = 2; % # of layers for optimal evaluation (To compensate the loss of solutions by weighted sum)

Db = zeros(Ndb,Nx);
db_ = RandDec( Nx,Nd );

D = zeros(0,Nx); EVs = []; EC = []; Nsol = zeros(Ndb,1);

for ii = 1:Ndb
    Db(ii,:) = db_;
    
    % Evaluate proxy measures for all Val(Dn), all Dn
    [EVs_prox, EC_, EVs_db] = EvalProxObj(S,S_,Px,Pm,Pl,Cd,db_);
    
    % Optimize in proxy space
    Dopt_ = Optimize( EVs_prox,EC_,Nlay,D );
    
    % Evaluate the obtained solution in original space
    if rem(ii,2)
        [EVs_, EC_, EVs_prox, AG] = EvalObj( Dopt_,S,S_,Px,Pm,Pl,Cd,EVs_prox,EVs_db,db_ );
    else
        [EVs_, EC_,EVs_prox] = EvalObj( Dopt_,S,S_,Px,Pm,Pl,Cd,EVs_prox,EVs_db,db_ );
    end
    subplot(1,2,1)
    semilogy( EVs_prox,EVs_,'sq' ); grid on; drawnow;
    
    % Sort non-dominated solutions in original problem
    [D,EVs,EC] = SortNonDominSol( [D; Dopt_],[EVs; EVs_],[EC;EC_] );
    Nsol(ii) = size(D,1);
    
    % Select Next db
    if rem(ii,2)
        db_ = selectNextBasis( EVs_,EVs_prox,Dopt_ );
    else
        db_ = RandDec( Nx,Nd );
    end      
    
    % Plot Optimal solutions
    subplot(1,2,2)
    semilogy( EC,EVs,'*' ); grid on; drawnow;  
    
    if ii == 2
        Fsz = 18; Fsz_label = 14; Msz = 9; Msz_diff = 2.5;

        Ngen = 1000; % # of generations of GA
        Niter = length( Nsol ); % # of iterations of the proposed method
        eval( ['load GA_SF_gen' num2str(Ngen)] )
        figure;
        semilogx( EVsga,ECga,'o','Markersize',Msz )
        hold on
        semilogx( EVs,EC,'*','Markersize',Msz )
        grid on
        axis([3e-4 4e-2 3.0e3 1e4])
%         ax = gca;
        ax.XAxis.FontSize = Fsz_label;
        ax.YAxis.FontSize = Fsz_label;
        str_leg_ga = ['Genetic Algorithm (' num2str(Ngen) ' gen.)'];
        str_leg_pro = [ 'Proposed method (' num2str(ii) ' iter.)' ];
        ylabel( '\Sigma {\itE} [{\it{V_n}} | {\it{d_n}} ]','Interpreter','tex','Fontsize',Fsz,'FontName','times new roman')
        xlabel('{\itE} [ {\it{V_{N+1}}} | {\bf{d}} ]','Interpreter','tex','Fontsize',Fsz,'FontName','times new roman')
        legend( {str_leg_ga str_leg_pro},...
            'Fontsize',Fsz-2,'FontName','times new roman' )

        str_fname_emf = [ 'figure/SF_result_' num2str(Ngen) 'gens.emf' ];
        str_fname_pdf = [ 'figure/SF_result_' num2str(Ngen) 'gens.pdf' ];
        saveas(gcf,str_fname_emf)
        saveas(gcf,str_fname_pdf)
        
        figure;
        semilogx( EVs_,EVs_prox,'*','Markersize',Msz,'MarkerEdgeColor',[0.8500 0.3250 0.0980] ); grid on; drawnow;
        ax = gca;
        ax.XAxis.FontSize = Fsz_label;
        ax.YAxis.FontSize = Fsz_label;
        xlabel( '{\itE} [ {\it{V_{N+1}}} | {\bf{d}} ]','Interpreter','tex','Fontsize',Fsz,'FontName','times new roman')
        ylabel('${\tilde{\it{E}}} [ \it{V_{N+1}} | {\bf{d}} ]$', 'Interpreter','latex','Fontsize',Fsz,'FontName','times new roman')
        
        str_leg2 = ['Solutions obtained at ' num2str(ii) 'nd iter.'];
        legend( {str_leg2},'Location','northwest','Fontsize',Fsz-2,'FontName','times new roman' )
                
        str_fname_emf2 = [ 'figure/SF_ex_prox_plot' num2str(ii) 'iter.emf' ];
        str_fname_pdf2 = [ 'figure/SF_ex_prox_plot' num2str(ii) 'iter.pdf' ];
        
        
        saveas(gcf,str_fname_emf2)
        saveas(gcf,str_fname_pdf2)

        eval(['save Optim_SF_' num2str(ii) 'iter D EVs EC Nsol Db'])
    end
    
    disp( [num2str(ii) '-th iter done.'] )
end

save Optim_SF D EVs EC Nsol Db

%% Evaluate upper bound of P(s0)
load Quant_SF2
[EVs_up, ~] = EvalObjGA2( D,S2,S2_,Path2,Px,Pm,Pl,Cd );
save('Optim_SF','EVs_up','-append')