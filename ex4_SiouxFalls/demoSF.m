% run Quant_SF
% run Optim_SF
% run Optim_SF2
% run GA_SF

load Quant_SF
load Optim_SF
load OPtim_SF2

%% Figure
Fsz = 18; Fsz_label = 14; Msz = 9; Msz_diff = 2.5;

Ngen = 1e4; % # of generations of GA
Niter = length( Nsol ); % # of iterations of the proposed method
eval( ['load GA_SF_gen' num2str(Ngen)] )
figure;
semilogx( EVsga,ECga,'o','Markersize',Msz )
hold on
semilogx( EVs,EC,'*','Markersize',Msz )
grid on
axis([3e-4 4e-2 3.0e3 1e4])
ax = gca;
ax.XAxis.FontSize = Fsz_label;
ax.YAxis.FontSize = Fsz_label;
str_leg_ga = ['Genetic Algorithm (' num2str(Ngen) ' gen.)'];
str_leg_pro = [ 'Proposed method (' num2str(Niter) ' iter.)' ];
xlabel( '{\itE} [ {\itV_{N+1}} | {\bf{\itd}} ]','Interpreter','tex','Fontsize',Fsz,'FontName','times new roman')
ylabel('\Sigma {\itE} [ {\itV_n} | {\itd_n} ]','Interpreter','tex','Fontsize',Fsz,'FontName','times new roman')
legend( {str_leg_ga str_leg_pro},...
    'Fontsize',Fsz-2,'FontName','times new roman' )

str_fname_emf = [ 'figure/SF_result_' num2str(Ngen) 'gens.emf' ];
str_fname_pdf = [ 'figure/SF_result_' num2str(Ngen) 'gens.pdf' ];
saveas(gcf,str_fname_emf)
saveas(gcf,str_fname_pdf)

%% Figure 2: Compare the results that use cut- and link-sets
[EVs_s,tmp] = sort( EVs ); EC_s = EC(tmp);
[EVs_up_s,tmp] = sort( EVs_up ); EC_up_s = EC(tmp);
[EVs2_s,tmp] = sort( EVs2 ); EC2_s = EC2(tmp);
[EVs2_low_s,tmp] = sort( EVs2_low ); EC_low_s = EC2(tmp);

line_w = 0.8; Msz2 = Msz-5.5;
figure;
semilogx( EVs_s,EC_s,'*--','Linewidth',line_w-.7,'Markersize',Msz2+2.1,'Color',[0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980] )
ax = gca;
ax.XAxis.FontSize = Fsz_label;
ax.YAxis.FontSize = Fsz_label;
axis([3e-4 4e-2 3.0e3 1e4])
hold on
grid on
semilogx( EVs2_low_s,EC_low_s,'o--','Linewidth',line_w,'Markersize',Msz2-.5,'Color',[0 0.4470 0.7410],'MarkerFaceColor',[0 0.4470 0.7410] )
semilogx( EVs_up_s,EC_up_s,'*--','Linewidth',line_w-.7,'Markersize',Msz2+2.1,'Color',[0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980] )
semilogx( EVs2_s,EC2_s,'o--','Linewidth',line_w,'Markersize',Msz2-.5,'Color',[0 0.4470 0.7410],'MarkerFaceColor',[0 0.4470 0.7410] )


xlabel( '{\itE} [ {\itV_{N+1}} | {\bf{\itd}} ]','Interpreter','tex','Fontsize',Fsz,'FontName','times new roman')
ylabel('\Sigma {\itE} [ {\itV_n} | {\itd_n} ]','Interpreter','tex','Fontsize',Fsz,'FontName','times new roman')
legend( {'Optimization of lower bound' 'Optimization of upper bound'},...
    'Fontsize',Fsz-2,'FontName','times new roman' )
saveas( gcf,'figure/SF_result_lowup.emf' )
saveas( gcf,'figure/SF_result_lowup.pdf' )