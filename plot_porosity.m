close all
clear all
load porosity590.mat
porosity=porosity/100;
index=depth<350;
depth=depth(index);
porosity=porosity(index);
[fitresult, gof]= createFit1_590(depth,porosity);

figure;
hold on
scatter(porosity,depth,1);
set(gca,'Ydir','reverse');
% h = plot( fitresult, depth, porosity );
f_b = fitresult.b;
f_a = 0.65-fitresult.b
f_c = fitresult.c;
save fp590.mat;
depth=sort(depth);
plot(f_a+f_b*exp(-depth/f_c),depth,'linewidth',2);

f_a=f_a-0.02;
plot(f_a+f_b*exp(-depth/f_c),depth,'--','linewidth',2);
f_a=f_a+0.04;
plot(f_a+f_b*exp(-depth/f_c),depth,'--','linewidth',2);

legend('measured porosity','exponential fitting','-0.02','+0.02','fontsize',8,'location','southeast');
% f_a = 0.525;
% f_b = 0.125;
% f_c = 175;
% plot(f_a+f_b*exp(-depth/f_c),depth,'linewidth',2);
xlabel('Porosity');
ylabel('Depth (m)');
title('Site 590');
set(gca,'FontSize',12,'FontWeight','bold')
% legend('measured porosity','our fitting','Richter and DePaolo (1987)','fontsize',8,'location','southeast');
print('porosity590.jpg','-djpeg','-r300');


clear all
[Leg,Site,Hi,Cor,T,Sc,Topcm,Depthmbsf,Densitygcc] = importfile_porosity('803.txt');
index=Densitygcc<2.7;
density=Densitygcc(index);
depth=Depthmbsf(index);
index=density>1;
density=density(index);
depth=depth(index);
porosity=(2.7-density)/1.7;
[fitresult, gof]= createFit(depth,porosity);

figure;
hold on
scatter(porosity,depth,1);
set(gca,'Ydir','reverse');
% h = plot( fitresult, depth, porosity );
f_a = fitresult.a;
f_b = fitresult.b;
f_c = fitresult.c;
save fp803.mat;
depth=sort(depth);
plot(f_a+f_b*exp(-depth/f_c),depth,'linewidth',2);

f_a=f_a-0.02;
plot(f_a+f_b*exp(-depth/f_c),depth,'--','linewidth',2);
f_a=f_a+0.04;
plot(f_a+f_b*exp(-depth/f_c),depth,'--','linewidth',2);

legend('measured porosity','exponential fitting','-0.02','+0.02','fontsize',8,'location','southeast');

% f_a = 0;
% f_b = 68.6/100;
% f_c = 1/0.00075;
% plot(f_a+f_b*exp(-depth/f_c),depth,'linewidth',2);
xlabel('Porosity');
ylabel('Depth (m)');
title('Site 803');
set(gca,'FontSize',12,'FontWeight','bold')
% legend('measured porosity','our fitting','Bassinot et al. (1993)','fontsize',8,'location','northwest');
print('porosity803.jpg','-djpeg','-r300');


clear all
[Leg,Site,Hi,Cor,T,Sc,Topcm,Depthmbsf,Densitygcc] = importfile_porosity('805.txt');
index=Densitygcc<2.7;
density=Densitygcc(index);
depth=Depthmbsf(index);
index=density>1;
density=density(index);
depth=depth(index);
porosity=(2.7-density)/1.7;
[fitresult, gof]= createFit(depth,porosity);

figure;
hold on
scatter(porosity,depth,1);
set(gca,'Ydir','reverse');
% h = plot( fitresult, depth, porosity );
f_a = fitresult.a;
f_b = fitresult.b;
f_c = fitresult.c;
save fp805.mat;
depth=sort(depth);
plot(f_a+f_b*exp(-depth/f_c),depth,'linewidth',2);

f_a=f_a-0.02;
plot(f_a+f_b*exp(-depth/f_c),depth,'--','linewidth',2);
f_a=f_a+0.04;
plot(f_a+f_b*exp(-depth/f_c),depth,'--','linewidth',2);

legend('measured porosity','exponential fitting','-0.02','+0.02','fontsize',8,'location','southeast');

% f_a = 0;
% f_b = 70.2/100;
% f_c = 1/0.00076;
% plot(f_a+f_b*exp(-depth/f_c),depth,'linewidth',2);
xlabel('Porosity');
ylabel('Depth (m)');
title('Site 805');
set(gca,'FontSize',12,'FontWeight','bold')
% legend('measured porosity','our fitting','Bassinot et al. (1993)','fontsize',8,'location','northwest');
print('porosity805.jpg','-djpeg','-r300');

clear all
[Leg,Site,Hi,Cor,T,Sc,Topcm,Depthmbsf,Densitygcc] = importfile_porosity('806.txt');
index=Densitygcc<2.7;
density=Densitygcc(index);
depth=Depthmbsf(index);
index=density>1;
density=density(index);
depth=depth(index);
porosity=(2.7-density)/1.7;
[fitresult, gof]= createFit(depth,porosity);

figure;
hold on
scatter(porosity,depth,1);
set(gca,'Ydir','reverse');
% h = plot( fitresult, depth, porosity );
f_a = fitresult.a;
f_b = fitresult.b;
f_c = fitresult.c;
save fp806.mat
depth=sort(depth);
plot(f_a+f_b*exp(-depth/f_c),depth,'linewidth',2);

f_a=f_a-0.02;
plot(f_a+f_b*exp(-depth/f_c),depth,'--','linewidth',2);
f_a=f_a+0.04;
plot(f_a+f_b*exp(-depth/f_c),depth,'--','linewidth',2);

legend('measured porosity','exponential fitting','-0.02','+0.02','fontsize',8,'location','southeast');

% f_a = 0;
% f_b = 69.2/100;
% f_c = 1/0.00045;
% plot(f_a+f_b*exp(-depth/f_c),depth,'linewidth',2);
xlabel('Porosity');
ylabel('Depth (m)');
title('Site 806');
set(gca,'FontSize',12,'FontWeight','bold')
% legend('measured porosity','our fitting','Bassinot et al. (1993)','fontsize',8,'location','northwest');
print('porosity806.jpg','-djpeg','-r300');



clear all
[Leg,Site,Hi,Cor,T,Sc,Topcm,Depthmbsf,Densitygcc] = importfile_porosity('807.txt');
index=Densitygcc<2.7;
density=Densitygcc(index);
depth=Depthmbsf(index);
index=density>1;
density=density(index);
depth=depth(index);
index=depth<400;
density=density(index);
depth=depth(index);
porosity=(2.7-density)/1.7;
[fitresult, gof]= createFit(depth,porosity);

figure;
hold on
scatter(porosity,depth,1);
set(gca,'Ydir','reverse');
% h = plot( fitresult, depth, porosity );
f_a = fitresult.a;
f_b = fitresult.b;
f_c = fitresult.c;
save fp807.mat
depth=sort(depth);
plot(f_a+f_b*exp(-depth/f_c),depth,'linewidth',2);

f_a=f_a-0.02;
plot(f_a+f_b*exp(-depth/f_c),depth,'--','linewidth',2);
f_a=f_a+0.04;
plot(f_a+f_b*exp(-depth/f_c),depth,'--','linewidth',2);

legend('measured porosity','exponential fitting','-0.02','+0.02','fontsize',8,'location','southeast');

% f_a = 0;
% f_b = 69.6/100;
% f_c = 1/0.00051;
% plot(f_a+f_b*exp(-depth/f_c),depth,'linewidth',2);
xlabel('Porosity');
ylabel('Depth (m)');
title('Site 807');
set(gca,'FontSize',12,'FontWeight','bold')
% legend('measured porosity','our fitting','Bassinot et al. (1993)','fontsize',8,'location','northwest');
print('porosity807.jpg','-djpeg','-r300');

