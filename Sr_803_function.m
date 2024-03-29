function C = Sr_803_function(zz,zzdata)

load('parameters_su_803.mat');
load('parameters_Ca_803.mat');

alpha=zz(1);
beta=zz(2);
gamma=zz(3);
v=zz(4);
gra_sr=zz(5);

%% ============define parameters==============================
K_sr_sea = 0.1935;    % Sr/Ca ratio in the solid over seawater

rho_s = 2.7;    % g/cm^3 ->g/m^3
rho_f = 1;      % g/cm^3 ->g/m^3

gra_su = x803(6);
v_su = x803(4);
gra_ca = y803(6);
v_ca = y803(5);

V=0;
Ks=x803(5);

load('workspace_803.mat');

%------- initialize the paramter------
J = length(z);
N = J;

% ----------Copy and Paste---------
fc = 89.30/100 ; 
% Porosity & Depth (fit required from initial report) 
f_a = 0;
f_b = 68.6/100;
f_c = 1/0.00075;

temp_T = 0.007926*z+275.1;
% ---------------------------------

phi_z = f_a + f_b*exp(-z/f_c);  % porosity over depth

R_d = zeros(N,J);
R_net = zeros(N,J);
R_p = zeros(N,J);

for n = 1:N
    for j = 1:n
        R_net(n,j) = 1*(y803(2)+y803(3)*exp(-z(n-j+1)/y803(4)));
    end
end

K_sr = zeros(length(z),1);

for i = 1:length(z)
    %     K_sr(i) = exp(-4.38 + 1.71e3/temp_T(i)-4.22e5/temp_T(i)^2)-0.006;
    K_sr(i)=0.02;
end

D0 = zeros(length(z),1);
D = zeros(length(z),1);

for i = 1:length(z)
    D0(i) = (3.69+ 0.169 * (temp_T(i)-273.15)) * 3.65 * 24 * 36;
    D(i) = D0(i) / (1-log(phi_z(i)^2));
end

L = 0.5;
k_su = x803(3);
G0 = x803(2);  %[mM/1 porewater]  by fitting

% K_sp is from fitting data, fitting celestite
K_sp = zeros(J,1);
for i = 1 :J
    K_sp(i) = -0.003557 * z(i) + 20.73+0.5;
end


for n = 1:N
    for j = 1:n
        R_d (n,j) = alpha+beta*exp(-z(n-j+1)/gamma);
        R_p(n,j) = R_d(n,j) - R_net(n,j);
    end
end

%% ==========Ca calculation and best fit===========

Ca_pore = zeros(N,J);
K_aq = 0;   % adsoprtion coefficient


for n =1:N
    Ca_pore(n,n:n+1) = 10.62 + 0.161*t(N-n+1);
end

% Ca_pore(:,1) = 10;
Ca_pore(2,1) = Ca_pore(2,2)+ z_step(2)* gra_ca;


%-----calculation
for n = 2: N-1
    
    a10 = 2*D(n-2+1)*t_step/z_step(n-2+2)^2;
    b10 = t_step*v_ca/z_step(n-2+2);
    
    Ca_pore(n+1,2) = ((a10-b10)*Ca_pore(n,2+1) + (a10+b10)*(Ca_pore(n, 3) + 2 * z_step(n)* gra_ca) + (1-a10)*Ca_pore(n-1,2) + 2*10*1000 * t_step * R_net(n-2+1) *fc* rho_s / rho_f * (1-phi_z(n-2+2)) / phi_z(n-2+2)) / (1+a10);
    
    for j = 3:n
        
        a1 = 2*D(n-j+1)*t_step/z_step(n-j+2)^2;
        b1 = t_step*v_ca/z_step(n-j+2);
        
        Ca_pore(n+1,j) = ((a1-b1)*Ca_pore(n,j+1) + (a1+b1)*Ca_pore(n,j-1) + (1-a1)*Ca_pore(n-1,j) + 2*10*1000 * t_step * R_net(n-j+1) *fc* rho_s / rho_f * (1-phi_z(n-j+2)) / phi_z(n-j+2)) / (1+a1);
        
    end
    
    Ca_pore(n+1, 1) = Ca_pore(n+1, 3) + 2 * z_step(n+1)* gra_ca;

end

% Ca_pore(:,1)=Ca_pore(:,2);



k = zeros(N,J);     % Sr/Sr2+ coefficient between solid and pore water

k_eq = 0.025;
K_sr = zeros(N,J);  % Sr/Ca coefficient
for n = 1: N
    for j = 1:n
        K_sr(n,j) = 0.24 / (1+ R_d(n-j+1)/R_p(n-j+1)*(0.24/k_eq-1));
    end
end

% K_sr(:,:) = 0.02;

for n = 1:N
    for j = 1: n
        k(n,j) = K_sr(n,j) * 10 / Ca_pore(n,j) * 1000;      % Sr/Ca ratio -> Sr ratio
    end
end

%% =========== Sulfate reduction model ==========

G= zeros(N,J);
G(:,1) = G0;
G(1,2) = G0;

for i = 1:N
    G(i,i:i+1) = G0;
end

C_su=zeros(N,J);

%---------boundary condition of sulfate in the pore water---------
Sea_sulfate=zeros(N,1);
for i = 1: N
    Sea_sulfate(i) = -0.15 * t(i) + 28;
end

for i = 1:N
    C_su(i,i:i+1) = Sea_sulfate(N-i+1);      % initial value of sulfate in the pore water
end

C_su(:,1) = Sea_sulfate(end);

%% ============ Sr and Sulfate calculation ===========

C_s = zeros(N, J); % in solid (ppm)
C_f = zeros(N, J); % in pore fluid (ppm)

% -----------Initial / Boundary Conditions-----------------

seawater803 = ones(1,N);
solid0 = zeros(1,N);

K = zeros(N,1);     % Sr/Sr2+ coefficient between solid and seawater
Ca_sea = zeros(N,1);

for n = 1: N
    Ca_sea(n) = 10.62+0.161*t(n);
end

for n = 1:N
    K(n) = K_sr_sea * 10 / Ca_sea(n)*1000;
end

% assume the initial condition of seawater
load lear_curve.mat
seawater803=transpose(lear_curve(t)).*Ca_sea/1000;

for i = 1:N
    solid0(i) = K (i) * seawater803(i);
end

for i = 1:N
    C_s(i,i:i+1) = solid0(N-i+1);  %initial value at the lastest deposit carbonate solid
end
C_s(1,2) = solid0(N);

for i = 1:N
    C_f(i, i:i+1) = seawater803(N-i+1);  %initial value at the lastest pore water
end
% C_f(2,1) = seawater803(N-1);
% C_f(:,1) = seawater803(N);
C_f(2,1) = C_f(2,2)+z_step(2)* gra_sr;

% -----Calculate the first concentration in solid and pore water-

for n = 2:N-1
    
    a10 = 2*D(n-2+1)*t_step/z_step(n-2+2)^2;
    b10 = t_step*v/z_step(n-2+2);
    b10_su = t_step*v_su/z_step(n-2+2);
    b20 = 2*V*t_step/z_step(n-2+2);
    
    a20 = 2*t_step*R_p(n-1,2)*rho_s/rho_f*(1-phi_z(n-2+2))/phi_z(n-2+2);  % attention the z(2) always represent the bottom
    a30 = R_d(n-1,2)*2*t_step;
    
    C_f(n,2-1) = C_f(n,2+1)+ 2*z_step(n)* gra_sr;
    C_su(n, 2-1) = C_su(n, 2+1)+ 2 * z_step(n)* gra_su;
    
    C_f(n+1,2) = ((a10-b10)*C_f(n,2+1)+(a10+b10)*C_f(n,2-1) + (1-a10-a20*fc*k(n-1,2))*C_f(n-1,2) + 2*t_step*R_d(n-1,2)*fc*rho_s/rho_f*(1-phi_z(n-2+2))/phi_z(n-2+2)*C_s(n-1,2)) / (1+a10);

    C_su(n+1,2) = ((a10-b10_su)*C_su(n,2+1)+(a10+b10_su)*C_su(n,2-1)+(1-a10)*C_su(n-1,2)-2*t_step*k_su*L*G(n-1,2)*C_su(n-1,2)/(C_su(n-1,2)+Ks))/(1+a10);
    G(n+1,2) = G(n-1,2) * (1 - 2*k_su*t_step*C_su(n-1,2)/(Ks +C_su(n-1,2))) - b20 * G(n,2+1) + b20 * G(n,2-1);
    
    if C_su(n+1,2) < 0
        C_su(n+1,2) = 0;
    end
    
    for j = 3:n
        
        a1 = 2*D(n-j+1)*t_step/z_step(n-j+2)^2;
        b1 = t_step*v/z_step(n-j+2);
        b1_su = t_step*v_su/z_step(n-j+2);
        b2 = 2*V*t_step/z_step(n-j+2);
        
        a2 = 2*t_step*R_p(n-1,j)*rho_s/rho_f*(1-phi_z(n-j+2))/phi_z(n-j+2);  % attention the z(2) always represent the bottom
        a3 = R_d(n-1,j)*2*t_step;
              
        C_f(n+1,j) = ((a1-b1)*C_f(n,j+1)+(a1+b1)*C_f(n,j-1) + (1-a1-a2*fc*k(n-1,j))*C_f(n-1,j) + 2*t_step*R_d(n-1,j)*fc*rho_s/rho_f*(1-phi_z(n-j+2))/phi_z(n-j+2)*C_s(n-1,j)) / (1+a1);
        C_s(n+1,j) = (1-a3)*C_s(n-1,j) + R_p(n-1,j)*2*t_step*k(n-1,j)*C_f(n-1,j)+b2*(C_s(n,j-1)-C_s(n,j+1));
        
        C_su(n+1,j) = ((a1-b1_su)*C_su(n,j+1)+(a1+b1_su)*C_su(n,j-1)+(1-a1)*C_su(n-1,j)-2*t_step*k_su*L*G(n-1,j)*C_su(n-1,j)/(C_su(n-1,j)+Ks))/(1+a1);
        G(n+1,j) = G(n-1,j) * (1 - 2*k_su*t_step*C_su(n-1,j)/(Ks +C_su(n-1,j))) - b2 * G(n,j+1) + b2 * G(n,j-1);
        
        if C_su(n+1,j) < 0
            C_su(n+1,j) = 0;
        end
%         % judge the precipitation of SrSO4
%         
%         if C_su(n+1,j) * C_f(n+1,j)  >=  K_sp(n+1-j+1)
%             C_f(n+1,j) = K_sp(n+1-j+1)/C_su(n+1,j);
%         end
    end
    
    C_f(n+1,1) = C_f(n+1,3) + 2*z_step(n+1)* gra_sr;
    C_s(n+1,1) = C_s(n+1,3);
    C_su(n+1,1) = C_su(n+1,3) + 2*z_step(n+1)* gra_su;
    
end
C=interp1(z(1:N),fliplr(C_f(N,1:J)),zzdata);

