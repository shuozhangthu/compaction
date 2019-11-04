% up to 6/22
% simulate the calcium in the pore water
close all;
% clear all;

Ca_807_main;

load('workspace_807.mat');

alpha=y(1);
beta=y(2);
gamma=y(3);
R_net = alpha+beta*exp(-z/gamma) ;

v=y(4);
gra_ca = y(5);

%------- initialize the paramter------
J = length(z);
N = J;

% ----------Copy and Paste---------
fc = 91.40/100 ; 

% Porosity & Depth (fit required from initial report) 
temp_T = 0.01633*z+274.9;

f_a = 0;
f_b = 69.6/100;
f_c = 1/0.00051;
phi_z = f_a + f_b*exp(-z/f_c);  % porosity over depth
% ---------------------------------

rho_s = 2.7;    % g/cm^3 ->g/m^3
rho_f = 1;      % g/cm^3 ->g/m^3

phi_z = f_a + f_b*exp(-z/f_c);  % porosity over depth

D0 = zeros(length(z),1);
D_ca = zeros(length(z),1);

for i = 1:length(z)
    D0(i) = (3.69+ 0.169 * (temp_T(i)-273.15)) * 3.65 * 24 * 36;
    D_ca(i) = D0(i) / (1-log(phi_z(i)^2));
end


%------- initialize the paramter------
J = length(z);
N = length(t);
Ca_pore = zeros(N,J);
K_aq = 0;   % adsoprtion coefficient 
% R_net = 10e-5 ;% R_net = R_d - R_p (dissoltion - precipitation)

% %% - Hide-
% k_r = 0.20e-5;
% 
% Ca_eq = zeros(N,J);
% 
% for n =1:N
%     Ca_eq(n,n) = 10.62 + 0.161 * t(N-n+1);
%     for j =1:n-1
% %         Ca_eq(n,j) = Ca_eq(n,n) + 0.05718 * z(n-j+1);
%         Ca_eq(n,j) = Ca_eq(n,n) - 2.976e-07 * z(n-j+1)^3 + 0.0002505 * z(n-j+1)^2 + 0.004831 * z(n-j+1);
%     end
% end
        
% R_net = 0.00035*ones(N,1);
% for i = 1:round(N/3)
%     R_net(i) = -0.0002;
% end


%---initial the Ca value, could paste to my code------
% gra_ca = 0.0118 ;

for n =1:N
    Ca_pore(n,n:n+1) = 10.62 + 0.161*t(N-n+1);    
end

% Ca_pore(:,1) = 10;
Ca_pore(2,1) = Ca_pore(2,2)+ z_step(2)* gra_ca;


%-----calculation
for n = 2: N-1
    
    a10 = 2*D_ca(n-2+1)*t_step/z_step(n-2+2)^2;
    b10 = t_step*v/z_step(n-2+2);

    Ca_pore(n+1,2) = ((a10-b10)*Ca_pore(n,2+1) + (a10+b10)*(Ca_pore(n, 3) + 2 * z_step(n)* gra_ca) + (1-a10)*Ca_pore(n-1,2) + 2*10*1000 * t_step * R_net(n-2+1) *fc* rho_s / rho_f * (1-phi_z(n-2+2)) / phi_z(n-2+2)) / (1+a10);
      
    for j = 3:n
        
        a1 = 2*D_ca(n-j+1)*t_step/z_step(n-j+2)^2;
        b1 = t_step*v/z_step(n-j+2);
        
        Ca_pore(n+1,j) = ((a1-b1)*Ca_pore(n,j+1) + (a1+b1)*Ca_pore(n,j-1) + (1-a1)*Ca_pore(n-1,j) + 2*10*1000 * t_step * R_net(n-j+1) *fc* rho_s / rho_f * (1-phi_z(n-j+2)) / phi_z(n-j+2)) / (1+a1);
        
    end
    
    Ca_pore(n+1, 1) = Ca_pore(n+1, 3) + 2 * z_step(n+1)* gra_ca;
%     Ca_pore(n, n+1) = Ca_pore(n+1, n);  % for the accuracy at the top ! no idea

%     Ca_pore(n, n+1) = Ca_pore(n+1, n+1);

%     Ca_pore(n, n+1) = Ca_pore(n, n);
    
end

% Ca_pore(:,1)=Ca_pore(:,2);

%% ----------Plot --------------

t_total = t(end);
%-------end of paste--------
t1 = round((t_total - 35)/t_step);
t2 = round((t_total - 25)/t_step);
t3 = round((t_total - 15)/t_step);

figure
plot(fliplr(Ca_pore(N,1:J)),z,'linewidth',2)
hold on
plot(fliplr(Ca_pore(N-t1,1:J-t1)),z(1+t1:end),'--','linewidth',2)
plot(fliplr(Ca_pore(N-t2,1:J-t2)),z(1+t2:end),':','linewidth',2)
plot(fliplr(Ca_pore(N-t3,1:J-t3)),z(1+t3:end),'-.','linewidth',2)
 
% 
% plot(fliplr(Ca_eq(N,:)),z)
% plot(fliplr(Ca_eq(N-165,1:J-165)),z(1+165:end))
% plot(fliplr(Ca_eq(N-282,1:J-282)),z(1+282:end))
% plot(fliplr(Ca_eq(N-399,1:J-399)),z(1+399:end))

set(gca,'Ydir','reverse')
title("Site 807");
xlabel('Calcium (mM) in pore water');
ylabel('Present-day depth (m)')
set(gca,'FontSize',12)

[Leg1,Site1,Topcm1,Botcm,Depthmbsf1,CalciumCamM,ChlorinityClmM,MagnesiumMgmM,pHpHna,SodiumNamM,StrontiumSruM,SulfateSO4mM,SilicaH4SiO4uM,AlkalinityALKmM,SalinitySALna] = importfile_water('water.xlsx');

index=(Site1==site_Number);
depth2=Depthmbsf1(index);
Ca2_data=CalciumCamM(index);

scatter (Ca2_data,depth2);

legend('t=46.3 m.y., present day','t=35 m.y.','t=25 m.y.','t=15 m.y.','Calcium Data')

print('calcium807.jpg','-djpeg','-r600');

% figure
% plot(R_ca,z,'--','linewidth',2);
% hold on
% plot(R_su_flip(2:end),z,'-','linewidth',2);
% set(gca,'Ydir','reverse')
% xlabel('R_{calcite} and R_{sulfate} (mM/Myr)');
% ylabel('Depth (m)')
% legend('R_{calcite}','R_{sulfate}');
% print('rates.jpg','-djpeg','-r600');

% %% save data
% 
% newName = 'z807';
% S.(newName) = z;
% 
% newName = 'Rca807';
% S.(newName) = R_ca;
% 
% newName = 'Rsu807';
% S.(newName) = R_su_flip(2:end);
% save('sulfate807.mat', '-struct', 'S'); 
