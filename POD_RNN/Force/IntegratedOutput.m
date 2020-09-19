% POF: 2020
% author: R.Gupta and S.R.Bukka 
% POD-RNN hybrid model 

%% this code carries integration over cylinder boundary and plots it
%% Note that use only pressure field here (not velocity)
% "Data_rigid_cylv3" containes finite-element mesh data 
% "P_recon.mat" is the reconstructed pressure field 
% "rigidV3.oisd" containes full-order pressure and viscous forces: 
% Forces from full order
% Total_x Total_y Pressure_x Pressure_y Viscous_x Viscous_y 

% make sure that you have reconstruvted pressure field "P_recon" generated
% from "reconstruction.m" in POD-RNN analysis

%%
clc
clear all
close all

%% load full-order finite-element data 
load ('Data_rigid_cylv3.mat');    % full-order mesh information: connectivity, coordinates, etc.   
load('../P_recon.mat');              % reconstructed pressure field 
cnn = conn;                       % connectivity 
crd = crd(:,2:4);                 % x,y,z field coordinates 
 
%Number of element nodes in 2D (4 for fourNodeQuad)
nen = 4 ;

%% Get FSI details 
% Swap BCCyl (it is 2 node element on cylinder surface)
% BCByl containes boundary nodes on the stationary cylinder 
temp = BCCyl(:,1);
BCCyl(:,1) = BCCyl(:,2);
BCCyl(:,2) = temp;
 
nElem = size(BCCyl,1);             % 156
ndof = size(unique(BCCyl),1) ;     % 156
 
% Get elements corresponding to first layer of boundary layer
[tf,idx2] = ismember(BCCyl(:,1),cnn(:,1));
[tf,idx2(:,2)] = ismember(BCCyl(:,1),cnn(:,2));
[tf,idx2(:,3)] = ismember(BCCyl(:,1),cnn(:,3));
[tf,idx2(:,4)] = ismember(BCCyl(:,1),cnn(:,4));
elemCyl = unique(idx2(:));
elemCyl(elemCyl==0)=[];
 
cnnCyl = cnn(elemCyl,:);
 
% Reorder the element points for reduced integration 
for i=1:size(elemCyl,1)
    [tf, idx1(i)] = ismember(BCCyl(i,1),cnnCyl(i,:));
    cnnCylNew(i,1) = cnnCyl(i,idx1(i));
     if (idx1(i) == 3)
         cnnCylNew(i,2) = cnnCyl(i,idx1(i)+1);
         cnnCylNew(i,3) = cnnCyl(i,1);
         cnnCylNew(i,4) = cnnCyl(i,2);
     elseif (idx1(i) == 4)
        cnnCylNew(i,2) = cnnCyl(i,1);
        cnnCylNew(i,3) = cnnCyl(i,2);
        cnnCylNew(i,4) = cnnCyl(i,3);
    end
end

%% do reduced numerical integration of Cauchy stress-tensor over cylinder boundary 

% Shape functions, gauss points and weights for Numerical integration
% Gauss points
gP = [-1/sqrt(3), -1.0
       1/sqrt(3), -1.0] ;
% Gauss weights
gW = [1.0, 1.0];
% Shape functions
N(:,1) = 0.25.*(1-gP(:,1)).*(1-gP(:,2)) ;   % just 2 here and rest is 0. 
N(:,2) = 0.25.*(1+gP(:,1)).*(1-gP(:,2)) ;
N(:,3) = 0.25.*(1+gP(:,1)).*(1+gP(:,2)) ;
N(:,4) = 0.25.*(1-gP(:,1)).*(1+gP(:,2)) ;
     
% Derivative of shape functions
Nx(:,1) = -0.25.*(1-gP(:,2)) ;
Nx(:,2) =  0.25.*(1-gP(:,2)) ;
Nx(:,3) =  0.25.*(1+gP(:,2)) ;
Nx(:,4) = -0.25.*(1+gP(:,2)) ;
Ny(:,1) = -0.25.*(1-gP(:,1)) ;
Ny(:,2) = -0.25.*(1+gP(:,1)) ;
Ny(:,3) =  0.25.*(1+gP(:,1)) ;
Ny(:,4) =  0.25.*(1-gP(:,1)) ;
Nx = Nx' ;
Ny = Ny' ;
% Number of quadrature points
nQuad = size(gW,2);

xxf = zeros(size(cnnCylNew));
yyf = zeros(size(cnnCylNew));

Force = zeros(2,1000); 

for k=1:1000                             % parent loop over all 1000 predicted time-steps
Sol.p = P_recon(:,k);     
    
% Localize the data to each element
 for i=1:nen
    xxf(:,i) = crd(cnnCylNew(:,i),1);
    yyf(:,i) = crd(cnnCylNew(:,i),2);
    pres(:,i) = Sol.p(cnnCylNew(:,i),1) ;
 end

for p = 1:nQuad  
    % Jacobian for each element
    J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
         yyf*[Nx(:,p)], yyf*[Ny(:,p)]];
    if size(J,2)==1
        J = J';
    end
    
    volume = sqrt(J(:,1).^2 + J(:,3).^2);
    vol = ( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
    vol = abs(vol);
    
    % Normal evaluation
    normal(:,1) = -J(:,3)./volume ;
    normal(:,2) = -J(:,1)./volume ;
    
    DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(vol,1,nen);
    DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(vol,1,nen);
    
    % Pressure 
    locP  = sum(repmat(N(p,:),nElem,1).*pres,2);
            
    % Length/ Area of the line/ surface integral in 1D/2D
    A0(:,p) = gW(p).*volume ;
    
    % X-traction 
    Fx_p(:,p) =  (-locP(:,1).*normal(:,1)).*A0(:,p);               
    
    % Y-traction                  
    Fy_p(:,p) = (-locP(:,1).*normal(:,2)) .*A0(:,p);
         
end

% Summation of all quadrature data
A0 = sum(A0,2);
Length = [sum(A0,1)];
Fx_p = sum(sum(Fx_p,2));
Fy_p = sum(sum(Fy_p,2));

% Write the force components here
Force(1,k) = Fx_p; 
Force(2,k) = Fy_p; 
clear Fx_p Fy_p 
end

%% Plot the calculated force 
normalization = 500.0; 
% Forces from full order 
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

% prediction data 
Force = Force'/normalization; 
Fx_pr = Force(:,1);  Fy_pr = Force(:,2); 

% Timing Details 
t_tr = [1:3500];            % training 
t_pr = [1:1000]';           % prediction  
t_fo_pr = [3501:4500]';     % full-order-prediction 
NA_x = 875:25:975; 

% lift force  
figure(1)
fo = Force_fo(3501:3900,4); pr = Fy_pr(1:400,1);  % in-phase 
plot((3501:3900)*0.25, pr,'r:','LineWidth',3.3); hold on; 
plot((3501:3900)*0.25, fo,'k-','LineWidth',2); 
xlabel('$\frac{tU_{\infty}}{D}$','Interpreter','latex');
ylabel('$\bf{C_{L,p}}$','Interpreter','latex');
set(gca,'linewidth',2.5); set(gca,'fontsize',32); 
set(gca, 'FontName', 'Times New Roman');
l = legend('prediction','actual'); hold on;
hold on; set(l, 'interpreter', 'latex'); axis square; 
title('Prediction','interpreter', 'latex'); 
set(gca, 'XTick',NA_x); hold on; 
axis([875 975 -0.8 0.8]);

% drag force
figure(2)
fo = Force_fo(3501:3900,3); pr = Fx_pr(1:400,1);  % in-phase 
plot((3501:3900)*0.25, pr,'r:','LineWidth',3.3); hold on; 
plot((3501:3900)*0.25, fo,'k-','LineWidth',2); 
xlabel('$\frac{tU_{\infty}}{D}$','Interpreter','latex');
ylabel('$\bf{C_{D,p}}$','Interpreter','latex');
set(gca,'linewidth',2.5); set(gca,'fontsize',32); 
set(gca, 'FontName', 'Times New Roman');
l = legend('prediction','actual'); hold on;
hold on; set(l, 'interpreter', 'latex'); axis square; 
title('Prediction','interpreter', 'latex'); 
set(gca, 'XTick',NA_x); hold on; 
axis([875 975 1.0 1.04]);