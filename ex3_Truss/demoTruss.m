% run Quant_Truss
% run Optim_Truss
% run GA_Truss

load Quant_Truss
load Optim_Truss

%% Figure
Fsz = 18; Fsz_label = 14; Msz = 9; Msz_diff = 2.5;

Ngen = 500; % # of generations of GA
Niter = length( Nsol ); % # of iterations of the proposed method
eval( ['load GA_Truss_gen' num2str(Ngen)] )
figure;
semilogx( EVsga,ECga,'o','Markersize',Msz )
hold on
semilogx( EVs,EC,'*','Markersize',Msz )
grid on
axis([1e-6 6e-1 0.55 1.15])
ax = gca;
ax.XAxis.FontSize = Fsz_label;
ax.YAxis.FontSize = Fsz_label;
str_leg_ga = ['Genetic Algorithm (' num2str(Ngen) ' gen.)'];
str_leg_pro = [ 'Proposed method (' num2str(Niter) ' iter.)' ];
ylabel( '\Sigma {\itE} [ {\itV_{n}} | {\itd_n} ]','Interpreter','tex','Fontsize',Fsz,'FontName','times new roman')
xlabel('{\itE} [ {\itV_{N+1}} | {\bf{d}} ]','Interpreter','tex','Fontsize',Fsz,'FontName','times new roman')
legend( {str_leg_ga str_leg_pro},...
    'Fontsize',Fsz-2,'FontName','times new roman' )

str_fname_emf = [ 'figure/Truss_result_' num2str(Ngen) 'gens.emf' ];
str_fname_pdf = [ 'figure/Truss_result_' num2str(Ngen) 'gens.pdf' ];
saveas(gcf,str_fname_emf)
saveas(gcf,str_fname_pdf)


