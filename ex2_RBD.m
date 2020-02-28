clear; close all;
import mbn.*; import func.*
rng(1)

%%
[Ncomp, Cpm, var, B, val, EVcost] = MBNquant_RBD;
[Pxd,CpmProd] = evalPxd( Cpm,var,B,var.S );
Nlay = 2;
Niter = 10;

%% Optimize
Ndec = cellfun( @(x) length(x.C),Cpm(var.D) );
basisDecRule = [];
decOpt = []; decOpt_EVsys = []; decOpt_EVcost = [];
for ii = 1:Niter
    
   if rem(ii,2)
       basisDecRule_i = getRandDec( Ndec );
   else
       [basisDecRule_i,basisDecRuleId_i] = selectNextBasis( decOpt_EVsys_EVcost_i,decOpt_prox_EVprox_EVcost_i,decOpt_prox_i );
       plot( decOpt_EVsys_EVcost_i(:,1),decOpt_prox_EVprox_EVcost_i(:,1),'o' ); hold on;
       plot( decOpt_EVsys_EVcost_i(basisDecRuleId_i,1),decOpt_prox_EVprox_EVcost_i(basisDecRuleId_i,1),'*' ); hold off; drawnow;
   end   
   basisDecRule = [basisDecRule;basisDecRule_i];
      
   proxEVsys_i = getProxEVSys( Cpm,var,basisDecRule_i,Pxd,CpmProd );
   EVcost_basis_i = evalEVcost( basisDecRule_i,EVcost );
   criticalWeight_i = getCriticalWeight( proxEVsys_i,EVcost );
   [decOpt_prox_i,decOpt_prox_EVprox_EVcost_i] = optimize( proxEVsys_i,EVcost,criticalWeight_i,Nlay );
   decOpt_EVsys_EVcost_i = evalEVsys_EVcost( decOpt_prox_i,CpmProd.p,Pxd,EVcost );
   [decOpt,decOpt_EVsys,decOpt_EVcost] = SortNonDominSol( [decOpt;decOpt_prox_i],[decOpt_EVsys; decOpt_EVsys_EVcost_i(:,1)],[decOpt_EVcost; decOpt_EVsys_EVcost_i(:,2)] );

end


%% All decision rules
Ndec = cellfun( @(x) length(x.C),Cpm(var.D) );
decAll = (1:Ndec(1))';
for nn = 2:Ncomp
    dec_n = (1:Ndec(nn))';
    decAll = [repmat( decAll,Ndec(nn),1 ) repelem( dec_n,size(decAll,1),1 )];
end
decAll_EVsys_EVcost = evalEVsys_EVcost( decAll,CpmProd.p,Pxd,EVcost );
[decAll_optim,decAll_optim_EVsys,decAll_optim_EVcost] = SortNonDominSol( decAll,decAll_EVsys_EVcost(:,1),decAll_EVsys_EVcost(:,2) );

%% figure
Fsz = 16; Fsz_label = 12; Msz = 9;
figure;
plot(decAll_EVsys_EVcost(:,1),decAll_EVsys_EVcost(:,2),'.')
hold on
plot(decAll_optim_EVsys,decAll_optim_EVcost,'sq','Markersize',Msz)
plot(decOpt_EVsys,decOpt_EVcost,'*','Markersize',Msz)
grid on

axis([0 0.023 650 2150])
ax = gca;
ax.XAxis.FontSize = Fsz_label;
ax.YAxis.FontSize = Fsz_label;
str_leg_all = 'All decision rules';
str_leg_ex = 'Non-dominated rules (exact)';
str_leg_pro = [ 'Proposed method (' num2str(Niter) ' iter.)' ];
xlabel( '{\itE} [ {\itV_{N+1}} | {\bf{d}} ]','Interpreter','tex','Fontsize',Fsz,'FontName','times new roman')
ylabel('\Sigma {\itE} [ {\itV_n} | {\it{d_n}} ]','Interpreter','tex','Fontsize',Fsz,'FontName','times new roman')
legend( {str_leg_all str_leg_ex str_leg_pro},...
    'Fontsize',Fsz-2,'FontName','times new roman' )

saveas(gcf,'figure/RBD_result.emf')
saveas(gcf,'figure/RBD_result.pdf')