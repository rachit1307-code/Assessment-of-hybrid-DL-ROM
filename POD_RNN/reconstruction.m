% POF: 2020
% author: R.Gupta and S.R.Bukka 
% POD-RNN hybrid model 

%% make sure you have run the script "pod.m" 
% either generate your own predictions using "lstmrecurrent.m" OR
% use the prediction that is already provided for reference in mat file: 
% "y_pred_n_pres.mat" : time-series prediction (pressure) 
% "y_pred_n_velx.mat" : time-series prediction (x-velocity)

%% Prediction
% load prediction 
load('y_pred_n_pres.mat'); 

y_pred_n = y_pred_n_pres; 
P_recon = zeros(26114,1000);
y_pred_n = y_pred_n';

% reconstruct the flow field 
for i = 1:1000
    P_recon(:,i) = Pmean' + Pmodes(:,:)*y_pred_n(:,i);
end

% get truth for comparion 
P_truth = zeros(26114,4500);
for i = 1:4500
    P_truth(:,i) = Xp(i,:)'+ Pmean';
end

P_error_local = zeros(26114,1);
P_error_rms_space = P_truth(:,3501:4500)-P_recon(:,:);

% spatial error 
for i = 1:1000
        P_error_recon(:,i) = abs(P_truth(:,3500+i)-P_recon(:,i))./norm(P_truth(:,3500+i),2);
end

% "P_recon" containes 1000 predicted fields from ts:[4001:5000]
% Reconstructed/predicted pressure field at t
figure(6)
cp_recon = griddata(X',Y',P_recon(:,800),Xc,Yc,'cubic');
contourf(Xc,Yc,cp_recon,10); hold on; axis square; colorbar; % caxis([-0.04 0.03]);
hold on; 
set(gca,'fontsize',32);
%colormap jet
hold on
fill(x,y,'w'); hold on; 
xlabel('X/D','Interpreter','Latex','FontSize',32); 
ylabel('Y/D','Interpreter','Latex','FontSize',32); 
title('Prediction','Interpreter','Latex','FontSize',32); hold on;


% True Pressure field at t 
figure(7)
P_truth_ct = griddata(X',Y',P_truth(:,4300),Xc,Yc,'cubic');
contourf(Xc,Yc,P_truth_ct,10); hold on; axis square; colorbar; % caxis([-0.04 0.03]);
hold on; 
set(gca,'fontsize',32);
%colormap jet
hold on
fill(x,y,'w'); hold on; 
xlabel('X/D','Interpreter','Latex','FontSize',32); 
ylabel('Y/D','Interpreter','Latex','FontSize',32); 
title('Truth','Interpreter','Latex','FontSize',32); hold on;

% error plot 
figure(8)
cp_recon4 = griddata(X',Y',P_error_recon(:,800),Xc,Yc);
contourf(Xc,Yc,cp_recon4,10); hold on; axis square; colorbar; caxis([0 5e-4]);
hold on; 
set(gca,'fontsize',32);
fill(x,y,'w');
xlabel('X/D','Interpreter','Latex','FontSize',32); 
ylabel('Y/D','Interpreter','Latex','FontSize',32); 
title('$E_{n}$','Interpreter','Latex','FontSize',32); hold on;

save('P_recon.mat','P_recon'); 