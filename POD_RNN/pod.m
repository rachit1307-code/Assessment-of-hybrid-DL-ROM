% POF: 2020
% author: R.Gupta and S.R.Bukka 
% POD-RNN hybrid model 

% this code generates the POD spatial and temporal modes from full-order
% data. 

clc
clear all 
close all 

%% load data 
% total snapshots = 5000
% snapshots for analysis = [501:5000]
% traning data = [501:3500]
% testing data = [3501:4500]
% ts = 0.25 s
ns = 4500;                               % time-snapshots 
nodes = 26114;                           % point cloud      
load('full_order_cyl_data/crd_rb.mat');  % coordinates 
load('full_order_cyl_data/pres_rb.mat'); % pressure  
load('full_order_cyl_data/vel_rb.mat');  % velocities  
X = (crd_rb(:,1))'; Y = (crd_rb(:,2))';
P = zeros(nodes,ns);                     % the field   


% select the snapshot data for analysis
P(:,:) = pres_rb(:,1,501:5000);      % select pressure or x-velocity data 
%P(:,:) = vel_rb(:,1,501:5000);

% mean and fluctuating matrix 
Pmean = zeros(nodes,1); 
for i=1:ns
    Pmean = Pmean + P(:,i); 
end    
Pmean = Pmean/ns; 
Xp = P - Pmean; 
clear vel_rb pres_rb crd_rb 
Pmean = Pmean'; Xp = Xp'; 

%% POD analysis
Nmodes = 5;         % number of modes
Rp_re = Xp*Xp';     % co-variance matrix
 
% eigen value of co-variance matric 
[Psip_re,Lp_re] = eig(Rp_re);   
Lp_re = abs(Lp_re);

% get total energy 
sumLp = 0;
for i=ns:-1:1
     sumLp = sumLp+Lp_re(i,i);     
end

% the modes with more than 95 percent energy!
sumLr = 0;
for i=ns:-1:1
     sumLr = sumLr+Lp_re(i,i);
     if sumLr/sumLp >0.95
         r=i;
         break
     end
end

% project and get the POD-modes 
Phip_re = (Xp'*Psip_re)*Lp_re^(-1/2); 
Pmodes = Phip_re(:,ns:-1:ns-Nmodes+1);

% select reference grid to post-process results 
delX = 10/63; delY = 10/63;   
X1 = [-5: delX :5]'; Y1 = [5: -delY: -5]'; 
[Xc,Yc] = meshgrid(X1,Y1); 

% do the cumulative energy analysis 
modes = [1:50]; 
CE = zeros(50,1);  CE(1,1) = Lp_re(ns,ns); 
k = 1; 
for i=(ns-1):-1:(ns-49)
    CE(k+1,1) = CE(k,1) + Lp_re(i,i);
    k = k+1;
end    
figure(1)
CE = CE/sumLp; 
h1 = plot(modes, CE*100, '-ro','LineWidth',3 ,'MarkerSize',12); hold on; 
h2 = plot(modes, 95*ones(50,1), '--k','LineWidth',3); 
axis square ;  
xlabel('Mode Number','Interpreter','latex');
ylabel('Cumulative Energy','Interpreter','latex');
title('Pressure','Interpreter','latex'); 
legend('Cumulative Energy','$95\%$ of total energy' ,'Interpreter','latex');
set(gca,'linewidth',2.5); set(gca,'fontsize',32); 
set(gca, 'FontName', 'Times New Roman'); hold on; 
axis([0 50 60 100]);

% do the total energy analysis 
TE = zeros(50,1); 
k = 1; 
for i=ns:-1:(ns-49)
    TE(k,1) = Lp_re(i,i);
    k = k+1; 
end     
figure(2)
TE = TE/sumLp;  
semilogy(modes, TE*100, '-bo','LineWidth',3 ,'MarkerSize',12); hold on; 
axis square; 
xlabel('Mode Number','Interpreter','latex');
ylabel('$\%$ Total Energy','Interpreter','latex');
title('Pressure','Interpreter','latex'); 
set(gca,'linewidth',2.5); set(gca,'fontsize',32); 
set(gca, 'FontName', 'Times New Roman'); hold on; 
axis([0 50 1.0e-15 1.0e+3]); 


figure(3)
mode = 1;                   % select spatial mode number to visualise 
cp = griddata(X,Y,Pmodes(:,mode),Xc,Yc, 'cubic'); 
contourf(Xc,Yc,cp,10); hold on; axis square; colorbar; 
hold on; 
set(gca,'fontsize',32);
hold on;
t = linspace(0,2*pi,100);   % mask out cylinder's surface 
x = 0.5*cos(t);
y = 0.5*sin(t);
fill(x,y,'w'); hold on; 
xlabel('X/D','Interpreter','Latex','FontSize',32); 
ylabel('Y/D','Interpreter','Latex','FontSize',32); 
title('Pressure: Mode 1','Interpreter','Latex','FontSize',32); hold on;


%% get the time coefficients!
a1 = zeros(Nmodes,ns);
for i = 1:Nmodes
    for j = 1:ns
         a1(i,j) = Xp(j,:)*(Pmodes(:,i));
     end
end

% plot how the tempral modes look like 
figure(4)

% time-history of mode 1
subplot(5,1,1)
plot((501:1500)*0.25, a1(1,1:1000),'k-','Linewidth',1.5); hold on; 
set(gca,'fontsize',16);
ylabel('$a_{1}$','Interpreter','latex','fontsize',32);
xlim([500*0.25 1500*0.25]);
title('Pressure','Interpreter','latex','fontsize',32);
set(gca, 'FontName', 'Times New Roman'); hold on;

% time-history of mode 2
subplot(5,1,2)
plot((501:1500)*0.25, a1(2,1:1000),'k-','Linewidth',1.5); hold on;
set(gca,'fontsize',16);
ylabel('$a_{2}$','Interpreter','latex','fontsize',32);
xlim([500*0.25 1500*0.25]);
set(gca, 'FontName', 'Times New Roman');
ax = gca; 
ax.YAxis.Exponent = 3; 

% time-history of mode 3
subplot(5,1,3)
plot((501:1500)*0.25, a1(3,1:1000),'k-','Linewidth',1.5); hold on;
set(gca,'fontsize',16);
ylabel('$a_{3}$','Interpreter','latex','fontsize',32);
xlim([500*0.25 1500*0.25]);
set(gca, 'FontName', 'Times New Roman');
ax = gca; 
ax.YAxis.Exponent = 3; 

% time-history of mode 4
subplot(5,1,4)
plot((501:1500)*0.25, a1(4,1:1000),'k-','Linewidth',1.5); hold on;
set(gca,'fontsize',16);
ylabel('$a_{4}$','Interpreter','latex','fontsize',32);
xlim([500*0.25 1500*0.25]);
set(gca, 'FontName', 'Times New Roman');
ax = gca; 
ax.YAxis.Exponent = 3; 

% time-history of mode 5
subplot(5,1,5)
plot((501:1500)*0.25, a1(5,1:1000),'k-','Linewidth',1.5); hold on;
set(gca,'fontsize',16);
ylabel('$a_{5}$','Interpreter','latex','fontsize',32);
xlabel('$\frac{tU_{\infty}}{D}$','Interpreter','latex','fontsize',32);
xlim([500*0.25 1500*0.25]);
set(gca, 'FontName', 'Times New Roman');
ax = gca; 
ax.YAxis.Exponent = 2; 
