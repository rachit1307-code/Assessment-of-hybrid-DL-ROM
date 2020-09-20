% POF: 2020
% author: R.Gupta 
% This utility predicts force on FSI using the pressure data from the observer-corrector method

%% pres files (Pressure) 
% rigidV3.oisd : contains the full-ouder force on the cylinder boundary for
% first 5000 timesteps: 
% Total_x Total_y Pressure_x Pressure_y Viscous_x Viscous_y 

% pres_rb_train_Fx/Fy: refers to the discrete force calculated from the
% training pressure field in x/y direction 
% pres_pred_32_Fx/Fy: refers to the discrete force calculated from the
% predicted pressure field in x/y direction 
% (refer to "utils.py" python script for it's calculation)

clc
clear all
close all

normalization = 500.0; 
%% Forces from full order (directly from CFD calculations)
% Total_x Total_y Pressure_x Pressure_y Viscous_x Viscous_y 
N = 5000; 
F1 = fopen('rigidV3.oisd','r');
for i=1:4
    fgets(F1); 
end    
Force_fo(1,:) = fscanf(F1,'%e %e %e %e %e %e',[1 6]);   
for i=2:N
    for j=1:7
        fgets(F1);
    end
    Force_fo(i,:) = fscanf(F1,'%e %e %e %e %e %e',[1 6]); 
end
Force_fo = Force_fo/normalization; 

%% training data (501:3500 time steps in full order, remember)
Fx_tr = -hdf5read('pres_train_Fx.h5','dataset')/normalization;     % x-training 
Fy_tr = -hdf5read('pres_train_Fy.h5','dataset')/normalization;     % y-training  

%% prediction data (3501:4500 time steps in full order, remember)
Fx_pr = -hdf5read('pres_pred_32_Fx.h5','decoded')/normalization;   % x-prediction
Fy_pr = -hdf5read('pres_pred_32_Fy.h5','decoded')/normalization;   % y-prediction 


%% x-traction prediction 

fo = Force_fo(501:900,3); tr = Fx_tr(1:400,1);  % in-phase 

% discrete force in x 
figure(1)
plot((501:900)*0.25, tr,'g:','LineWidth',3); hold on; 
plot((501:900)*0.25,fo,'k-','LineWidth',2); 
xlabel('$\frac{tU_{\infty}}{D}$','Interpreter','latex');
ylabel('$\bf{C_{D,p}}$','Interpreter','latex');
set(gca,'linewidth',2.5); set(gca,'fontsize',28); 
ax = gca;
set(gca, 'FontName', 'Times New Roman');
l = legend('discrete-force','full-order'); hold on;
hold on; set(l, 'interpreter', 'latex'); axis square; 
title('Training','interpreter', 'latex'); 
axis([124 225 1 1.35]);

% observer-corrector algorithm 
midd_tr = mean(tr); 
tr = tr - midd_tr; 

midd_fo = mean(fo); 
fo = fo - midd_fo;

Ef = max(fo)/max(tr); 
                                        
tr = tr*Ef;                                      % mapping function

tr = tr + midd_fo; 
fo = fo + midd_fo ;

% mapped discrete force using calculated Ef
figure(2)
plot((501:900)*0.25, tr,'g:','LineWidth',3); hold on; 
plot((501:900)*0.25,fo,'k-','LineWidth',2); 
xlabel('$\frac{tU_{\infty}}{D}$','Interpreter','latex');
ylabel('$\bf{C_{D,p}}$','Interpreter','latex');
set(gca,'linewidth',2.5); set(gca,'fontsize',28); 
ax = gca;
set(gca, 'FontName', 'Times New Roman');
l = legend('mapped-discrete-force','full-order'); hold on;
hold on; set(l, 'interpreter', 'latex'); axis square; 
title('Training','interpreter', 'latex'); 
axis([124 225 1 1.04]);


% x-traction prediction using predicted discrete force and Ef
pr = Fx_pr(1:400,1);     
pr = pr - mean(pr); 
pr = pr*Ef + midd_fo; 

NA_x = 875:25:975; 

figure(3)
plot((3501:3900)*0.25, pr,'r:','LineWidth',3.3); hold on; 
plot((3501:3900)*0.25, Force_fo(3501:3900,3),'k-','LineWidth',2); 
xlabel('$\frac{tU_{\infty}}{D}$','Interpreter','latex');
ylabel('$\bf{C_{D,p}}$','Interpreter','latex');
set(gca,'linewidth',2.5); set(gca,'fontsize',32); 
ax = gca;
set(gca, 'FontName', 'Times New Roman');
l = legend('prediction','actual'); hold on;
hold on; set(l, 'interpreter', 'latex'); axis square; 
title('Prediction','interpreter', 'latex'); 
set(gca, 'XTick',NA_x); hold on; 
axis([875 975 1.0 1.04]);

%% y-traction prediction 
fo = Force_fo(501:900,4); tr = Fy_tr(1:400,1);  % in-phase 

% discrete force in y
figure(4)
plot((501:900)*0.25, tr,'g:','LineWidth',3); hold on; 
plot((501:900)*0.25,fo,'k-','LineWidth',2); 
xlabel('$\frac{tU_{\infty}}{D}$','Interpreter','latex');
ylabel('$\bf{C_{L,p}}$','Interpreter','latex');
set(gca,'linewidth',2.5); set(gca,'fontsize',28); 
ax = gca;
set(gca, 'FontName', 'Times New Roman');
l = legend('discrete-force','full-order'); hold on;
hold on; set(l, 'interpreter', 'latex'); axis square; 
title('Training','interpreter', 'latex'); 
axis([124 225 -0.8 0.8]);

% observer-corrector algorithm 
midd_tr = mean(tr); 
tr = tr - midd_tr; 

midd_fo = mean(fo); 
fo = fo - midd_fo;

Ef = max(fo)/max(tr); 
                                        
tr = tr*Ef;                                      % mapping function

tr = tr + midd_fo; 
fo = fo + midd_fo ;

% mapped discrete force using calculated 
figure(5)
plot((501:900)*0.25, tr,'g:','LineWidth',3); hold on; 
plot((501:900)*0.25,fo,'k-','LineWidth',2); 
xlabel('$\frac{tU_{\infty}}{D}$','Interpreter','latex');
ylabel('$\bf{C_{L,p}}$','Interpreter','latex');
set(gca,'linewidth',2.5); set(gca,'fontsize',28); 
ax = gca;
set(gca, 'FontName', 'Times New Roman');
l = legend('mapped-discrete-force','full-order'); hold on;
hold on; set(l, 'interpreter', 'latex'); axis square; 
title('Training','interpreter', 'latex'); 
axis([124 225 -0.8 0.8]);


% y-traction prediction!
pr = Fy_pr(1:400,1);     
pr = pr - mean(pr); 
pr = pr*Ef + midd_fo; 

NA_x = 875:25:975; 

figure(6)
plot((3501:3900)*0.25, pr,'r:','LineWidth',3.3); hold on; 
plot((3501:3900)*0.25, Force_fo(3501:3900,4),'k-','LineWidth',2); 
xlabel('$\frac{tU_{\infty}}{D}$','Interpreter','latex');
ylabel('$\bf{C_{L,p}}$','Interpreter','latex');
set(gca,'linewidth',2.5); set(gca,'fontsize',32); 
ax = gca;
set(gca, 'FontName', 'Times New Roman');
l = legend('prediction','actual'); hold on;
hold on; set(l, 'interpreter', 'latex'); axis square; 
title('Prediction','interpreter', 'latex'); 
set(gca, 'XTick',NA_x); hold on; 
axis([875 975 -0.8 0.8]);
