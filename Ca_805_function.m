function C=Ca_805_function(y,ydata)

alpha=y(1);
beta=y(2);
gamma=y(3);
v=y(4);
gra_ca = y(5);

load('workspace_805.mat');

%------- initialize the paramter------
J = length(z);
N = J;

% ----------Copy and Paste---------
fc = 91.50/100 ; 

% Porosity & Depth (fit required from initial report) 
temp_T = 0.01651*z+275.8;

f_a = 0;
f_b = 70.2/100;
f_c = 1/0.00076;
phi_z = f_a + f_b*exp(-z/f_c);  % porosity over depth

% ---------------------------------

rho_s = 2.7;    % g/cm^3 ->g/m^3
rho_f = 1;      % g/cm^3 ->g/m^3

R_net=alpha+beta*exp(-z/gamma);

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


C=interp1(z(1:N),fliplr(Ca_pore(N,1:J)),ydata);
