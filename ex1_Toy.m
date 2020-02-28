clear; close all;
import mbn.*; import func.*

%%
[Ncomp, Cpm, var, B, val, EVcost] = MBNquant_Toy;
[Pxd,CpmProd] = evalPxd( Cpm,var,B,[var.S var.H] );
Nlay = 1;

%% Optimize
basisDecRule = ones( 1,Ncomp );
proxEVsys = getProxEVSys( Cpm,var,basisDecRule,Pxd,CpmProd );
EVcost_basis = evalEVcost( basisDecRule,EVcost );
criticalWeight = getCriticalWeight( proxEVsys,EVcost );
[decOpt_prox,decOpt_prox_EVprox_EVcost] = optimize( proxEVsys,EVcost,criticalWeight,Nlay );
decOpt_EVsys_EVcost = evalEVsys_EVcost( decOpt_prox,CpmProd.p,Pxd,EVcost );

%% All decision rules
Ndec = cellfun( @(x) length(x.C),Cpm(var.D) );
decAll = (1:Ndec(1))';
for nn = 2:Ncomp
    dec_n = (1:Ndec(nn))';
    decAll = [repmat( decAll,Ndec(nn),1 ) repelem( dec_n,size(decAll,1),1 )];
end
decAll_EVsys_EVcost = evalEVsys_EVcost( decAll,CpmProd.p,Pxd,EVcost );
[decAll_optim,decAll_optim_EVsys,decAll_optim_EVcost] = SortNonDominSol( decAll,decAll_EVsys_EVcost(:,1),decAll_EVsys_EVcost(:,2) );

%% Deciding the next decision rule
D = -diff( decOpt_prox_EVprox_EVcost(:,1) ) ./ diff( decOpt_EVsys_EVcost(:,1) );
[~,nextBasisIdx] = max(D); 
disp(['Next basis rule: ' num2str( decOpt_prox(nextBasisIdx,:) )])

%% figure
Fsz = 16; Fsz_label = 12; Msz = 9; lineWidth = 1;
figure;
plot(decAll_optim_EVsys,decAll_optim_EVcost,'sq')
hold on
plot(decAll_optim_EVsys,decAll_optim_EVcost,'*')
plot(decAll_optim_EVsys(setdiff(1:5,3)),decAll_optim_EVcost(setdiff(1:5,3)),'--','Linewidth',lineWidth)
grid on

ax = gca;
ax.XAxis.FontSize = Fsz_label;
ax.YAxis.FontSize = Fsz_label;
xlabel('{\itE} [ {\itV_{N+1}} | {\bf{d}} ]','Interpreter','tex','Fontsize',Fsz,'FontName','times new roman')
ylabel('\Sigma {\itE} [ {\itV_n} | {\itd_n} ]','Interpreter','tex','Fontsize',Fsz,'FontName','times new roman')

legend( {'All decision rules' 'Non-dominated solutions' 'Solutions by weighted sum'},'Fontsize',Fsz,'FontName','times new roman' )

saveas(gcf,'figure/Ex_Comp3_result.emf')
saveas(gcf,'figure/Ex_Comp3_result.pdf')