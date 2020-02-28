function Pex = Ind_Fragility_Ac(Sa,PGVv,d_limit,Ac)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This fuction calculates an exeedance probability for a specific damage %
% state. Input parameters are: Sa: spectral acceleration, PGV: Peak     % 
% Ground Velocity, d_limit: drift limit(0.02 for above the insignificat %
% level, 0.04 for above the moderate level, default for above the heavy %
% level)                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import SF_fragility.*
%% Outputs
Pex = zeros(2,1);
    % Pex(1): Pristine Pex for above the given level
    % Pex(2): Pristine Pex for above the heavy level

%% Inputs
g=386.4;
gm_group = 1.0;
dc_det = log(d_limit);   % d_limit: drift limit
path_p = ['PushData_34_pristine'];
path_d = ['PushData_34_damaged'];

figure_show = 0.0;   % when figure_show = 1, then the plots are given
                     % when figure_show = 0, the valuse of fragilities is given
A_r = log([linspace(0.01,0.6,10) linspace(0.7,0.9,5) linspace(1.0,4,10)]);
PGV_r1 = log(linspace(0.001,3,20));
PGV_r2 = PGV_r1;

SA = Sa*g;     % sec/in^2 %% From the relationship between Sa and PGA(PGA, from the PGA attenuation relationship)!!!
PGV = PGVv/2.54;  % in/sec   %% From the PGV relationship!!! 

%% Misc.

pushfile = [path_p '/bilinear.mat'];  % obtain Vy, stiffness, disp_y & alpha
massfile = [path_p '/constrated_mass.out'];
load (pushfile)
load (massfile)
T1=2*pi*(constrated_mass/stiffness)^0.5;    % Natural Period

pushfile = [path_d '/bilinear.mat'];  % obtain Vy, stiffness, disp_y & alpha
massfile = [path_d '/constrated_mass.out'];
load (pushfile)
load (massfile)
T2=2*pi*(constrated_mass/stiffness)^0.5;    % Natural Period

x = load ('ExampleBridge_34.txt');
soil_num = fix(x(1,11));               % soil type
Hc = x(1,4)*12;                        % column height

if figure_show == 0
    A_r = log(SA/g);
    PGV_r1 = log(PGV*T1/Hc);
    PGV_r2 = log(PGV*T2/Hc);
end

m=length(A_r);
n=length(PGV_r1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% capacity %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% prinstine
mu_dc(1) = -2.8276;       % drift capacity [log(dh)+gamma] obtained from "section_drift.tcl"
mu_vc(1) = 0.0806;     % shear capacity [log(vh)+gamma] obtained from "section_shear.tcl"

sigma_dc = 0.415; 
sigma_vc = 0.196;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% demand %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_dd0 = -2.7963;
theta_dd1 = 0.2771;
theta_dd2 = -0.9472;
theta_dd3 = 0.4575;
theta_dd4 = -0.1498;
theta_dd5 = 0.5417;

theta_vd0 = -0.5510;
theta_vd1 = -0.1071;
theta_vd2 = -0.2655;
theta_vd3 = 0.1652;
theta_vd4 = -0.1941;
theta_vd5 = 0.1773;

sigma_dd0 = 0.4306; 
sigma_vd0 = 0.3559;

wi1_d = 1.0;
wi2_d = 1.0;
wi1_v = 1.120960;
wi2_v = 9.019452;

%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% bi-variant %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho1 = -0.463;
rho2 = 0.6226;
mu = [0 0];

%% fragility: 1st case
n1 = 0; n2 = 0;
PGV_r = PGV_r1;
for k=1:m;
    [d_hat1,shear_hat1,Hc,ftAc,T,disp_y,Vy,alpha] = GetDemand_N2_Ac(exp(A_r(k)),path_p, Ac);
    dh(k) = d_hat1/Hc;
    vh(k) = shear_hat1/ftAc;
    dy = disp_y/Hc;
    vy = Vy/ftAc;
    Sd_r(k)=exp(A_r(k))*386.4/(2*pi/T)^2/Hc;
    fctor(k) = (sign(log(dh(k)/dy))+1)/2;
        
    if dh(k)<dy
        n1 = n1 + 1;
    else
        n2 = n2 + 1;
    end
    
    for j=1:n;      
        
        mu_dd(k,j) = log(dh(k)) + theta_dd0 + theta_dd1*log(vh(k)) + theta_dd2*log(Sd_r(k)) ...
            +  theta_dd3*PGV_r(j) + theta_dd4*gm_group + theta_dd5*(log(dh(k))-log(dy))*fctor(k);

        mu_vd(k,j)= theta_vd0+(1+theta_vd1)*log(vh(k))+theta_vd2*log(Sd_r(k))+theta_vd3*PGV_r(j)...
            +theta_vd4*(log(dh(k))-log(dy))*(1-fctor(k))+theta_vd5*(log(dh(k))-log(dy))*fctor(k);  

        if dh(k)<dy
            sigma_dd = sigma_dd0/wi1_d^0.5;
            sigma_vd = sigma_vd0/wi1_v^0.5;
        else 
            sigma_dd = sigma_dd0/wi2_d^0.5;
            sigma_vd = sigma_vd0/wi2_v^0.5;
        end
        sigma_d=(sigma_dc^2+sigma_dd^2)^0.5;
        sigma_v=(sigma_vc^2+sigma_vd^2)^0.5;
        rho = (rho1*sigma_dc*sigma_vc+rho2*sigma_dd*sigma_vd)/sigma_d/sigma_v;
        sigma=[1 rho; rho 1];
        
        gd(k,j)=(mu_dc(1)-mu_dd(k,j))/sigma_d;
        fd(k,j)=1-normcdf(gd(k,j),0,1);       
        gv(k,j)=(mu_vc(1)-mu_vd(k,j))/sigma_v;    
        fv(k,j)=1-normcdf(gv(k,j),0,1);
        
        X=[gd(k,j) gv(k,j)];       
        fbi(k,j) = fd(k,j)+fv(k,j)-mvncdf(-X,mu,sigma);
        
        sigma_det = [1 0; 0 1];
        gd_det(k,j)=(dc_det-mu_dd(k,j))/sigma_d;
        fd_det(k,j)=1-normcdf(gd_det(k,j),0,1);
        X_det=[gd_det(k,j) gv(k,j)];       
        fbi_det(k,j) = fd_det(k,j)+fv(k,j)-mvncdf(-X_det,mu,sigma_det);     
        
        if dh(k)<dy
            fd1_case1(n1,j)=fd(k,j);            
            fv1_case1(n1,j)=fv(k,j);            
            fbi1_case1(n1,j)=fbi(k,j); 
            A_r1_case1(n1)=A_r(k);
            fbi_det1_case1(n1,j)=fbi_det(k,j);
        else
            fd2_case1(n2,j)=fd(k,j);             
            fv2_case1(n2,j)=fv(k,j);           
            fbi2_case1(n2,j)=fbi(k,j);
            A_r2_case1(n2)=A_r(k);
            fbi_det2_case1(n2,j)=fbi_det(k,j);
        end         
    end   
       
end

if figure_show == 0
    % disp('---- pristine: bi-variant fragility with ultimate capacity ----')
    % disp(fbi)
    Pex(2)=fbi;
    % disp('---- pristine: bi-variant fragility with given service capacity ----')
    % disp(fbi_det)
    Pex(1)=fbi_det;
end

%% figures

if figure_show == 1.0
    
    %% figures: probabilistic capacity and deterministic capacity (pristine)
    xup=3;
    yup=4;
    x=linspace(0.1,0.9,9);
    name_det = num2str(exp(dc_det)*100);

    %%%%%%% drift-shear %%%%%%%
    scrsz = get(0,'ScreenSize');
    figure('Position',[100 100 scrsz(4)/2.8 scrsz(4)/3.0])
    [C,h] = contour(exp(PGV_r),exp(A_r1_case1),fbi1_case1,x,'LineColor',[0 0 0],'LineStyle','-','linewidth', 1.0);
    hold on;
    [C,h] = contour(exp(PGV_r),exp(A_r1_case1),fbi_det1_case1,x,'LineColor',[0 0 0],'LineStyle',':','linewidth', 1.0);
    axis([0,xup,0,yup])
    legend('ultimate',[name_det '% ']);legend('boxoff');
    % clabel(C,h,'manual')
    hold on;

    [C,h] = contour(exp(PGV_r),exp(A_r2_case1),fbi2_case1,x,'LineColor',[0 0 0],'LineStyle','-','linewidth', 1.0);
    hold on;
    [C,h] = contour(exp(PGV_r),exp(A_r2_case1),fbi_det2_case1,x,'LineColor',[0 0 0],'LineStyle',':','linewidth', 1.0);
    axis([0,xup,0,yup])
    % clabel(C,h,'manual')
    ylabel('PSA/g','FontName','Times New Roman');
    xlabel('PGV\timesT_1/H_c','FontName','Times New Roman');
    axis([0,xup,0,yup])

    %% figures: probabilistic capacity and deterministic capacity (damaged)

    %%%%%%% drift-shear %%%%%%%
    scrsz = get(0,'ScreenSize');
    figure('Position',[100 100 scrsz(4)/2.8 scrsz(4)/3.0])
    [C,h] = contour(exp(PGV_r),exp(A_r1_case2),fbi1_case2,x,'LineColor',[0 0 0],'LineStyle','-','linewidth', 1.0);
    hold on;
    [C,h] = contour(exp(PGV_r),exp(A_r1_case2),fbi_det1_case2,x,'LineColor',[0 0 0],'LineStyle',':','linewidth', 1.0);
    axis([0,xup,0,yup])
    legend('ultimate',[name_det '% ']);legend('boxoff');
    % clabel(C,h,'manual')
    hold on;

    [C,h] = contour(exp(PGV_r),exp(A_r2_case2),fbi2_case2,x,'LineColor',[0 0 0],'LineStyle','-','linewidth', 1.0);
    hold on;
    [C,h] = contour(exp(PGV_r),exp(A_r2_case2),fbi_det2_case2,x,'LineColor',[0 0 0],'LineStyle',':','linewidth', 1.0);
    axis([0,xup,0,yup])
    % clabel(C,h,'manual')
    ylabel('PSA/g','FontName','Times New Roman');
    xlabel('PGV\timesT_1/H_c','FontName','Times New Roman');
    axis([0,xup,0,yup])
end

end
