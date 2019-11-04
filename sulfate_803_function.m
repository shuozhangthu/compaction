function C = sulfate_function_803(x,xdata)

G0 = x(1);  %[mM/1 porewater]  by fitting
k_su=x(2);
v=x(3);
Ks = x(4);
gra_su= x(5);

% G and concentration in surfate
load('workspace_803.mat');

L = 0.5;
J = length(z_step);
N = J;

% t_step = 1/w;

V=0;

C_su=zeros(N,J);

% ----------Copy and Paste---------
fc = 89.30/100 ; 
% Porosity & Depth (fit required from initial report) 
f_a = 0;
f_b = 68.6/100;
f_c = 1/0.00075;

temp_T = 0.007926*z+275.1;

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

C=interp1(z(1:N),fliplr(C_su(N,1:J)),xdata);




