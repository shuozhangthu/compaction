% G and concentration in surfate
clear all
close all

sulfate_805_main

G0 = x(1);  %[mM/1 porewater]  by fitting
k_su=x(2);
v=x(3);
Ks = x(4);
gra_su=x(5);

load('workspace_805.mat');

L = 0.5;
J = length(z_step);
N = J;

V=0;

C_su=zeros(N,J);

% ----------Copy and Paste---------
fc = 91.50/100 ; 
% Porosity & Depth (fit required from initial report) 
f_a = 0;
f_b = 70.2/100;
f_c = 1/0.00076;

temp_T = 0.01651*z+275.8;
% ---------------------------------

phi_z = f_a + f_b*exp(-z/f_c);  % porosity over depth

D0 = zeros(length(z),1);
D_ca = zeros(length(z),1);
for i = 1:length(z)
    D0(i) = (3.69+ 0.169 * (temp_T(i)-273.15)) * 3.65 * 24 * 36;
    D_ca(i) = D0(i) / (1-log(phi_z(i)^2));
end

% =========== boundary condition============

Sea_sulfate = zeros(N,1);

for i = 1: N
    Sea_sulfate(i) = -0.15 * t(i) + 28;
end

for i = 1:N
    C_su(i,i:(i+1)) = Sea_sulfate(N-i+1);      % initial value of sulfate in the pore water
end

C_su(:,1) = 0;

G=zeros(N,J);

G(:,1) = G0;
G(1,2) = G0;

for i = 1:N
    G(i,i:i+1) = G0;
end

D = 7500;     % cm^2/s ->m^2/million year

for n = 2:N-1
    
    a10 = 2*D_ca(n-2+1)*t_step/z_step(n-2+2)^2;
    
    b10 = t_step*v/z_step(n-2+2);
    b20 = 2*V*t_step/z_step(n-2+2);
    
    C_su(n+1,2) = ((a10-b10)*C_su(n,2+1)+(a10+b10)*(C_su(n, 3)+ 2 * z_step(n)* gra_su)+(1-a10)*C_su(n-1,2)-2*t_step*k_su*L*G(n-1,2)*C_su(n-1,2)/(C_su(n-1,2)+Ks))/(1+a10);
    
    G(n+1,2) = G(n-1,2) * (1 - 2*k_su*t_step*C_su(n-1,2)/(Ks +C_su(n-1,2)));
    
    if C_su(n+1,2) < 0
        C_su(n+1,2) = 0;
    end
    
    for j = 3:n
        
        a1 = 2*D_ca(n-j+1)*t_step/z_step(n-j+2)^2;
        
        b1 = t_step*v/z_step(n-j+2);
        b2 = 2*V*t_step/z_step(n-j+2);
        
        C_su(n+1,j) = ((a1-b1)*C_su(n,j+1)+(a1+b1)*C_su(n,j-1)+(1-a1)*C_su(n-1,j)-2*t_step*k_su*L*G(n-1,j)*C_su(n-1,j)/(C_su(n-1,j)+Ks))/(1+a1);
        
        G(n+1,j) = G(n-1,j) * (1 - 2*k_su*t_step*C_su(n-1,j)/(Ks +C_su(n-1,j))) - b2 * G(n,j+1) + b2 * G(n,j-1);
        
        if C_su(n+1,j) < 0
            C_su(n+1,j) = 0;
        end
        
    end
    
end
C_su(:,1)=C_su(:,2);
G(:,1)=G(:,2);

[Leg1,Site1,Topcm1,Botcm,Depthmbsf1,CalciumCamM,ChlorinityClmM,MagnesiumMgmM,pHpHna,SodiumNamM,StrontiumSruM,SulfateSO4mM,SilicaH4SiO4uM,AlkalinityALKmM,SalinitySALna] = importfile_water('water.xlsx');

t_total = t(end);
t1 = round((t_total - 20)/t_step);
t2 = round((t_total - 15)/t_step);
t3 = round((t_total - 10)/t_step);

figure;
plot(fliplr(C_su(N,1:J)),z(1:N),'linewidth',2)
hold on
plot(fliplr(C_su(N-t1,1:J-t1)),z(1+t1:end),'--','linewidth',2)
plot(fliplr(C_su(N-t2,1:J-t2)),z(1+t2:end),':','linewidth',2)
plot(fliplr(C_su(N-t3,1:J-t3)),z(1+t3:end),'-.','linewidth',2)
set(gca,'Ydir','reverse')
title("Site 805");
xlabel('Sulfate (mM) in pore water');
ylabel('Present-day depth (m)')
set(gca,'FontSize',12)

index3=(Site1==site_Number & SulfateSO4mM >0);
depth3=Depthmbsf1(index3);
sulfate_data=SulfateSO4mM(index3);

scatter(sulfate_data,depth3,'ko')
legend('t=27.9 m.y., present day','t=20 m.y.','t=15 m.y.','t=10 m.y.','Sulfate Data','Location','SouthEast')
print('sulfate805.jpg','-djpeg','-r600');

figure;
plot(fliplr(C_su(N,1:J)),z(1:N),'linewidth',2)
hold on
plot(fliplr(C_su(N-t1,1:J-t1)),z(1:N-t1),'--','linewidth',2)
plot(fliplr(C_su(N-t2,1:J-t2)),z(1:N-t2),':','linewidth',2)
plot(fliplr(C_su(N-t3,1:J-t3)),z(1:N-t3),'-.','linewidth',2)
set(gca,'Ydir','reverse')
title('Sulfate concentration in pore water vs depth')
xlabel('Sulfate (mM) in pore water');
ylabel('Depth (m)')
set(gca,'FontSize',12)

index3=(Site1==site_Number & SulfateSO4mM >0);
depth3=Depthmbsf1(index3);
sulfate_data=SulfateSO4mM(index3);

scatter(sulfate_data,depth3,'ko')
legend('t=27.9 m.y., present day','t=20 m.y.','t=15 m.y.','t=10 m.y.','Sulfate Data','Location','SouthEast')
% print('sulfate2.jpg','-djpeg','-r600');

figure;
plot(fliplr(G(N,1:J)),z(1:N),'linewidth',2)
hold on
plot(fliplr(G(N-t1,1:J-t1)),z(1+t1:end),'--','linewidth',2)
plot(fliplr(G(N-t2,1:J-t2)),z(1+t2:end),':','linewidth',2)
plot(fliplr(G(N-t3,1:J-t3)),z(1+t3:end),'-.','linewidth',2)
set(gca,'Ydir','reverse')
% title('Organic carbon vs depth')
xlabel('G (mM)');
ylabel('Depth (m)')
set(gca,'FontSize',12)

legend('t=27.9 m.y., present day','t=20 m.y.','t=15 m.y.','t=10 m.y.','Location','SouthEast')
print('G805.jpg','-djpeg','-r600');

