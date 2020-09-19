% POF: 2020
% author: R.Gupta and S.R.Bukka 
% POD-RNN hybrid model 

%% make sure you have run "pod.m" and havepredicted modes saved in mat file
% this code is just a pos-process to pot predicted vs actual modes 
% load the predctied temporal modes from "lstmrecurrent.m" OR 
% use the prediction already provided for reference: 
% "y_pred_n_pres.mat" : time-series prediction (pressure) 
% "y_pred_n_velx.mat" : time-series prediction (x-velocity)

%% plot the predicted modes 
load('y_pred_n_pres.mat'); 
pred = y_pred_n_pres; 

figure(5)
subplot(5,1,1)
plot((3501:4500)*0.25, pred(1:1000,1),'r:','Linewidth',3); hold on; 
plot((3501:4500)*0.25, a1(1,3001:4500),'k-','Linewidth',1.5); hold on; 
set(gca,'fontsize',16);
ylabel('$a_{1}$','Interpreter','latex','fontsize',32);
xlim([3500*0.25 4500*0.25]);
legend('Prediction','Truth','Interpreter','Latex','Orientation','horizontal','fontsize',32);
title('Pressure','Interpreter','latex','fontsize',32);
set(gca, 'FontName', 'Times New Roman'); hold on;


subplot(5,1,2)
plot((3501:4500)*0.25, pred(1:1000,2),'r:','Linewidth',3); hold on;
plot((3501:4500)*0.25, a1(2,3501:4500),'k-','Linewidth',1.5); hold on;
set(gca,'fontsize',16);
ylabel('$a_{2}$','Interpreter','latex','fontsize',32);
xlim([3500*0.25 4500*0.25]);
set(gca, 'FontName', 'Times New Roman');
ax = gca; 
ax.YAxis.Exponent = 3; 

subplot(5,1,3)
 plot((3501:4500)*0.25, pred(1:1000,3),'r:','Linewidth',3); hold on;
plot((3501:4500)*0.25, a1(3,3501:4500),'k-','Linewidth',1.5); hold on;
set(gca,'fontsize',16);
ylabel('$a_{3}$','Interpreter','latex','fontsize',32);
xlim([3500*0.25 4500*0.25]);
set(gca, 'FontName', 'Times New Roman');
ax = gca; 
ax.YAxis.Exponent = 3; 
 

subplot(5,1,4)
plot((3501:4500)*0.25, pred(1:1000,4),'r:','Linewidth',3); hold on;
plot((3501:4500)*0.25, a1(4,3501:4500),'k-','Linewidth',1.5); hold on;
set(gca,'fontsize',16);
ylabel('$a_{4}$','Interpreter','latex','fontsize',32);
xlim([3500*0.25 4500*0.25]);
set(gca, 'FontName', 'Times New Roman');
ax = gca; 
ax.YAxis.Exponent = 3; 
 

subplot(5,1,5)
plot((3501:4500)*0.25, pred(1:1000,5),'r:','Linewidth',3); hold on;
plot((3501:4500)*0.25, a1(5,3501:4500),'k-','Linewidth',1.5); hold on;
set(gca,'fontsize',16);
ylabel('$a_{5}$','Interpreter','latex','fontsize',32);
xlim([3500*0.25 4500*0.25]);
xlabel('$\frac{tU_{\infty}}{D}$','Interpreter','latex','fontsize',32);
set(gca, 'FontName', 'Times New Roman');
ax = gca; 
ax.YAxis.Exponent = 2; 
 

