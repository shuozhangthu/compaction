clc
clear all;
close all;

[Leg2,Site2,Depthmbsf2,AgeMa] = importfile_age('age.xlsx');
load age.mat;

figure;
hold on

index=(Site2==803);
age_file=AgeMa(index);
depth_file=Depthmbsf2(index);

scatter(age_file,depth_file,'o');

index=(Site2==805);
age_file=AgeMa(index);
depth_file=Depthmbsf2(index);

scatter(age_file,depth_file,'s');

index=(Site2==806);
age_file=AgeMa(index);
depth_file=Depthmbsf2(index);

scatter(age_file,depth_file,'p');

index=(Site2==807);
age_file=AgeMa(index);
depth_file=Depthmbsf2(index);

scatter(age_file,depth_file,'d');

scatter(ans(:,1),ans(:,2),'h');

set(gca,'Ydir','reverse')

xlabel('Age (My)');
ylabel('Present-day depth (m)')
set(gca,'FontSize',12)
legend('Site 803','805','806','807','590B','Location','NorthEast');
box on

print('age_depth.jpg','-djpeg','-r600');



load('data803.mat');
load('data805.mat');
load('data806.mat');
load('data807.mat');
% load('rittenhouse1.mat');
% load('rittenhouse2.mat');

t_803=0.0591;
t_805=0.0385;
t_806=0.0271;
t_807=0.0329;

z807=z807(1:707);
Re807=Re807(1:707);
Rd807=Rd807(1:707);

t_590=1/38;
z590=transpose([0:1:500]);
Rd590=transpose(0.0027+0.067*exp(-z590/38/2.1));
Re590 = cumsum(Rd590)*t_590*100; % in percentage (%)
Re803c = cumsum(Rd803)*t_803*100; % in percentage (%)
Re805c = cumsum(Rd805)*t_805*100; % in percentage (%)
Re806c = cumsum(Rd806)*t_806*100; % in percentage (%)
Re807c = cumsum(Rd807)*t_807*100; % in percentage (%)

% f_a = 0.525;
% f_b = 0.125;
% f_c = 175;
load fp590.mat;

phi590=f_a+f_b*exp(-z590/f_c);   % porosity corresponding to the depth grid
comp_rate_590=(phi590-f_a)./(1-phi590)/f_c/t_590;

% 
% f_a = 0;
% f_b = 68.6/100;
% f_c = 1/0.00075;
load fp803.mat;

phi803=f_a+f_b*exp(-z803/f_c);   % porosity corresponding to the depth grid
comp_rate_803=(phi803-f_a)./(1-phi803)/f_c/t_803;


% f_a = 0;
% f_b = 70.2/100;
% f_c = 1/0.00076;
load fp805.mat;

phi805=f_a+f_b*exp(-z805/f_c);   % porosity corresponding to the depth grid
comp_rate_805=(phi805-f_a)./(1-phi805)/f_c/t_805;

% f_a = 0;
% f_b = 69.2/100;
% f_c = 1/0.00045;
load fp806.mat;

phi806=f_a+f_b*exp(-z806/f_c);   % porosity corresponding to the depth grid
comp_rate_806=(phi806-f_a)./(1-phi806)/f_c/t_806;

% 
% f_a = 0;
% f_b = 69.6/100;
% f_c = 1/0.00051;
load fp807.mat;

phi807=f_a+f_b*exp(-z807/f_c);   % porosity corresponding to the depth grid
comp_rate_807=(phi807-f_a)./(1-phi807)/f_c/t_807;


Ce803 = cumsum(Rd803.*transpose(1-phi803))*t_803*100; % in percentage (%)
Ce805 = cumsum(Rd805.*transpose(1-phi805))*t_805*100; % in percentage (%)
Ce806 = cumsum(Rd806.*transpose(1-phi806))*t_806*100; % in percentage (%)
Ce807 = cumsum(Rd807.*transpose(1-phi807))*t_807*100; % in percentage (%)
Ce590 = cumsum(Rd590.*transpose(1-phi590))*t_590*100; % in percentage (%)


figure;
hold on

plot(phi803,z803,'-','linewidth',2);
plot(phi805,z805,'-.','linewidth',2);
plot(phi806,z806,':','linewidth',2);
plot(phi807,z807,'--','linewidth',2);
plot(phi590,z590,'-','linewidth',2);


% set(gca, 'XScale', 'log')
set(gca,'Ydir','reverse')

xlabel('Porosity');
ylabel('Present-day depth (m)')
set(gca,'FontSize',12)
legend('Site 803','805','806','807','590B','Location','SouthEast');
box on

print('porosity_depth.jpg','-djpeg','-r600');


figure
hold on

plot(Rd803*100,z803,'-','linewidth',2);
plot(Rd805*100,z805,'-.','linewidth',2);
plot(Rd806*100,z806,':','linewidth',2);
plot(Rd807*100,z807,'--','linewidth',2);
plot(Rd590*100,z590,'-','linewidth',2);

set(gca,'Ydir','reverse')

xlabel('Recrystallization rate (%/Myr)');
ylabel('Depth (m)')
set(gca,'FontSize',12)
legend('Site 803','805','806','807','590B','Location','SouthEast');
box on 

print('Rd_depth.jpg','-djpeg','-r600');





figure;
hold on

% plot(reaction_rate_803*100,z_803,'r-','linewidth',2);
% plot(reaction_rate_805*100,z_805,'b-','linewidth',2);
% plot(reaction_rate_806*100,z_806,'g-','linewidth',2);
% plot(reaction_rate_807*100,z_807,'k-','linewidth',2);

plot(comp_rate_803*100,z803,'-','linewidth',2);
plot(comp_rate_805*100,z805,'-.','linewidth',2);
plot(comp_rate_806*100,z806,':','linewidth',2);
plot(comp_rate_807*100,z807,'--','linewidth',2);
plot(comp_rate_590*100,z590,'-','linewidth',2);



% set(gca, 'XScale', 'log')
set(gca,'Ydir','reverse')

xlabel('Porosity reduction rate (vol%/Myr)');
ylabel('Present-day depth (m)')
set(gca,'FontSize',12)
legend('Site 803','805','806','807','590B','Location','SouthEast');
box on

print('porosity_reduction_rate_depth.jpg','-djpeg','-r600');




figure;
hold on

plot(Re803,z803,'-','linewidth',2);
plot(Re805,z805,'-.','linewidth',2);
plot(Re806,z806,':','linewidth',2);
plot(Re807,z807,'--','linewidth',2);
plot(Re590,z590,'-','linewidth',2);
% plot(Re803c,z803,'-','linewidth',2);
% plot(Re805c,z805,'-.','linewidth',2);
% plot(Re806c,z806,':','linewidth',2);
% plot(Re807c,z807,'--','linewidth',2);

h=annotation('arrow');
set(h,'parent', gca, 'position', [23 294 2 0]);

h=annotation('arrow');
set(h,'parent', gca, 'position', [12.5 217 2 0]);

h=annotation('arrow');
set(h,'parent', gca, 'position', [15 339 2 0]);

h=annotation('arrow');
set(h,'parent', gca, 'position', [12 293 2 0]);

h=annotation('arrow');
set(h,'parent', gca, 'position', [13.5 269 2 0]);
% set(gca, 'XScale', 'log')
set(gca,'Ydir','reverse')

xlabel('Fraction of recrystallized calcite (%)');
ylabel('Present-day depth (m)')
set(gca,'FontSize',12)
legend('Site 803','805','806','807','590B','Location','SouthEast');
box on

print('recrystallized.jpg','-djpeg','-r600');




figure;
hold on

scatter(20,Rd803(1)*100,40,'filled','o');
scatter(18,Rd805(1)*100,40,'filled','s');
scatter(25,Rd806(1)*100,80,'filled','p');
scatter(29,Rd807(1)*100,40,'filled','d');
scatter(10,Rd590(1)*100,40,'filled','h');

xlabel('Mean grain size (\mum)');
ylabel('Recrystallization rate (vol%/Myr)')
set(gca,'FontSize',12)
legend('Site 803','805','806','807','590B','Location','SouthEast');
box on

print('grain_recry.jpg','-djpeg','-r600');

figure;
hold on

scatter(comp_rate_803(1)*100,Rd803(1)*100,40,'filled','o');
scatter(comp_rate_805(1)*100,Rd805(1)*100,40,'filled','s');
scatter(comp_rate_806(1)*100,Rd806(1)*100,80,'filled','p');
scatter(comp_rate_807(1)*100,Rd807(1)*100,40,'filled','d');
scatter(comp_rate_590(1)*100,Rd590(1)*100,40,'filled','h');

x=[0:1:12];
plot(x,x);
% plot(x,0.8*x);
% plot(x,0.6*x);
h1=text(6,6,'R=C');
set(h1,'Rotation',45);
% 
% h2=text(3.2,2.7,'R=0.8C');
% set(h2,'Rotation',40);
% 
% h3=text(4,2.5,'R=0.6C');
% set(h3,'Rotation',35);
% 
xlim([0,12]);
ylim([0,12]);

xlabel('Porosity reduction rate (vol%/Myr)');
ylabel('Recrystallization rate (vol%/Myr)')
set(gca,'FontSize',12)
legend('Site 803','805','806','807','590B','Location','SouthEast');
box on

print('porosity_reduction_rate_recrystallization.jpg','-djpeg','-r600');
% 
% figure;
% hold on
% 
% plot(comp_rate_803*100,Rd803*100,'-','linewidth',2);
% plot(comp_rate_805*100,Rd805*100,'-.','linewidth',2);
% plot(comp_rate_806*100,Rd806*100,':','linewidth',2);
% plot(comp_rate_807*100,Rd807*100,'--','linewidth',2);
% plot(comp_rate_590*100,Rd590*100,'-','linewidth',2);
% % set(gca, 'XScale', 'log')
% % set(gca, 'YScale', 'log')
% 
% x=[0:0.1:8];
% plot(x,x,'k-.','linewidth',1);
% % plot(x,0.8*x);
% % plot(x,0.6*x);
% h1=text(4,4.2,'R=C');
% set(h1,'Rotation',45);
% 
% % h2=text(3.2,2.7,'R=0.8C');
% % set(h2,'Rotation',40);
% % 
% % h3=text(4,2.5,'R=0.6C');
% % set(h3,'Rotation',35);
% % % 
% xlim([0,8]);
% ylim([0,8]);
% 
% xlabel('Total compaction rate (vol%/Myr)');
% ylabel('Recrystallization rate (vol%/Myr)')
% set(gca,'FontSize',12)
% legend('Site 803','805','806','807','590B','Location','SouthEast');
% box on
% 
% print('compaction_recrystallization2.jpg','-djpeg','-r600');
% 
% 
% figure;
% hold on
% 
% plot(phi803(1)*100-phi803*100,Ce803,'-','linewidth',2);
% plot(phi805(1)*100-phi805*100,Ce805,'-.','linewidth',2);
% plot(phi806(1)*100-phi806*100,Ce806,':','linewidth',2);
% plot(phi807(1)*100-phi807*100,Ce807,'--','linewidth',2);
% plot(phi590(1)*100-phi590*100,Ce590,'-','linewidth',2);
% x=[0:1:30];
% y=x;
% plot(x,y,'k-.','linewidth',1);
% % 
% % plot(rittenhouse1(:,1)+rittenhouse1(:,2),rittenhouse1(:,2),'k--','linewidth',2);
% % plot(rittenhouse2(:,1)+rittenhouse2(:,2),rittenhouse2(:,2),'k--','linewidth',2);
% % scatter(7.64,1.1);
% xlabel('Total porosity loss (%)');
% ylabel('Cement (%)')
% set(gca,'FontSize',12)
% 
% legend('Site 803','805','806','807','590B','1:1 curve','Location','best');
% box on
% 
% print('cement_porosity_loss.jpg','-djpeg','-r600');