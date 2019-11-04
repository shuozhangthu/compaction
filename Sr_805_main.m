clc
clear all
close all

site_Number = 805;


[Leg1,Site1,Topcm1,Botcm,Depthmbsf1,CalciumCamM,ChlorinityClmM,MagnesiumMgmM,pHpHna,SodiumNamM,StrontiumSruM,SulfateSO4mM,SilicaH4SiO4uM,AlkalinityALKmM,SalinitySALna] = importfile_water('water.xlsx');
index=(Site1==site_Number & StrontiumSruM>0 &Depthmbsf1<160);
depth=Depthmbsf1(index);
Sr=StrontiumSruM(index);

[depth, a_order] = sort(depth);
Sr = Sr(a_order,:);
alpha0 = 0;  %[mM/1 porewater]  by fitting
beta0=0.013;
gamma0=600;
v0=0;
gra_sr0=0.008;

zz0=[alpha0,beta0,gamma0,v0,gra_sr0];

% lb=[0 ,-1,477.28,0,0];
% ub=[0,1,477.28,0,0];
lb=[0 ,-1,408.2,0,0];
ub=[0,1,408.2,0,0];

% zz= zz0
zz = lsqcurvefit(@Sr_805_function,zz0,depth,Sr/1000,lb,ub)

figure;
plot(Sr_805_function(zz,depth),depth,'linewidth',2)
hold on
scatter(Sr/1000,depth,'ko')
hold on
set(gca,'Ydir','reverse')
title("Sr concentration in pore water (site:805)")
xlabel('Sr (mM) in pore water');
ylabel('Depth (m)')
set(gca,'FontSize',12)

newName = 'zz805';
S.(newName) = [site_Number,zz];
save('parameters_sr_805.mat', '-struct', 'S'); 


newName = 'fit_sr_805';
S.(newName) = Sr_805_function(zz,depth);
save('fit_sr_805.mat', '-struct', 'S');