% clc;
% clear all;
% close all;
% % A_r = linspace(log(0.01),log(10),20);
% % A_r = exp(A_r);
% A_r = 1.0;

function [d_hat1,shear_hat1,Hc,ftAc,T,disp_y,Vy,alpha_hardening] = GetDemand_N2_Ac(A_r,path,Ac)

g=386.4; 
groupfile = 'ExampleBridge_34.txt';
x = load (groupfile);
soil_num = fix(x(1,11));               % soil type
Hc = x(1,4)*12;                        % column height
fc = x(1,7);
Ac = Ac*(39.370079)^2; % m2 -> in2
ft = 6*(fc*1000)^0.5/1000;
ftAc = ft*Ac;

% Sa = Sa_r * g;
% Sd = Sd_r * Hc;

for groupnum = 34:34;
 
    pushfile = [path '/bilinear.mat'];  % obtain Vy, stiffness, disp_y & alpha
    massfile = [path '/constrated_mass.out'];
    load (pushfile)
    load (massfile)
    
    m=constrated_mass;
    w=m*g;
    k=stiffness;
    T=2*pi*(m/k)^0.5;
    omega = (k/m)^0.5;
    beta_damping=0.05;

%% Calculate the displacement and shear using deterministic procedural model (N2)
    
    A=A_r*386.4;
    Sd=A/omega^2;
    Ay=Vy/m;
    Ry1=max(A/Ay,1);                     % yield rduction factor

    x_alpha = [0 0.02 0.10];
    y_a = [1.0 1.0 0.8];
    y_b = [0.42 0.37 0.29];

    a = interp1(x_alpha,y_a,alpha_hardening,'linear','extrap'); 
    b = interp1(x_alpha,y_b,alpha_hardening,'linear','extrap');
    c = T^a/(1+T^a) + b/T;
    mu = 1 + (Ry1.^c-1)./c;

    if Sd <= disp_y
        d_hat1 = Sd;
        shear_hat1 = d_hat1*stiffness;
    else        
        d_hat1 = mu*disp_y;
        shear_hat1 = (d_hat1-disp_y)*stiffness*alpha_hardening + Vy;
    end
    
%     shear_hat1 = interp1(displacement,force,d_hat1,'linear','extrap'); 
    
%     d_hat1/Hc
%     shear_hat1/ftAc
% figure(1)
% plot(displacement,force)
% figure(2)
% plot(A_r,shear_hat1)

end

