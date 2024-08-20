% aps

clc; close all; clear;
% B-field 


saveDir = checkDir('/Users/rishabhdatta/Dropbox (MIT)/PUFFIN/Data/MARZ/bdot/z3697/processed/');

fig = openfig([saveDir, 'bdot_voltage.fig'])
legend('5 mm (A)', '5 mm (B)', '8 mm', '11 mm (A)', '11 mm (B)', '15 mm');

saveas(gcf, [saveDir, 'bdot_voltage-aps.png']); 

fig = openfig([saveDir, 'bfield.fig'])
% legend('5 mm (1)', '5 mm (2)', '8 mm', '11 mm (1)', '11 mm (2)', '15 mm','Current');
legend('off');
yyaxis left
ylim([0,35])
yyaxis right
ylim([0,35])
ax = gca;
ax.YAxis(2).Color = [0.5,0.5,0.5];

saveas(gcf, [saveDir, 'bfield-aps.png']); 

 yyaxis left
 magnifyPlot([220-60,220+60],[0,5],[.18 .5 .4 .35],2,true,[0.3,0.3,0.3]);
 xticks(220-60:20:220+60)
 grid on;
     
saveas(gcf, [saveDir, 'bfield-aps-nsf.png']); 

fig = openfig([saveDir, 'bfield-aps-2.fig'])
% legend('5 mm (1)', '5 mm (2)', '8 mm', '11 mm (1)', '11 mm (2)', '15 mm','Current');
legend('off');

ax = gca;
ax.YAxis(2).Color = [0.5,0.5,0.5];

% data = csvread('/Users/rishabhdatta/Dropbox (MIT)/PUFFIN/Data/MARZ/bdot/B_field.csv'); 
% 
% yyaxis left
% plot(data(:,1),data(:,2),'--','color',[0 0.5 0],'Linewidth',3,'DisplayName','Simulation')

data = dlmread('/Users/rishabhdatta/Dropbox (MIT)/PUFFIN/Data/MARZ/bdot/rad_trans_sim.csv',' '); 

yyaxis left

tid = linspace(0,500,size(data,1));
plot(tid,data(:,end),'-.','color',[0 0.5 0],'Linewidth',4,'DisplayName','Simulation')

ylim([0,35])
set(gca, 'Fontsize',32);
% formatPlots(900,20)

yyaxis right
ylim([0,35])

saveas(gcf, [saveDir, 'bfield-aps-3.png']); 



%%
% velocity


fig = openfig([saveDir, 'transit_time_nsf.fig'])
legend('8 mm','11 mm');
xlim([0,600]);
ylim([-200,500]);
set(gcf,'Position',[0 0 900 450]);
yticks([-200,0,200,400])
saveas(gcf, [saveDir, 'trasit_time-nsf.png']); 



fig = openfig([saveDir, 'velocity_displacement_nsf.fig'])
legend('Average Velocity','Linear Fit');
xlim([0,600]);
saveas(gcf, [saveDir, 'velocity_displacement-aps.png']); 

% data = csvread('C:\Users\rdatta\Dropbox (MIT)\PUFFIN\Data\MARZ\bdot\velocity.csv'); 
% 
% plot(data(:,1),data(:,2)/1e3,'--r','Linewidth',3,'DisplayName','Simulation')
set(gcf,'Position',[0 0 900 450]);
 ylim([0,150]);
saveas(gcf, [saveDir, 'velocity_displacement-2-nsf.png']); 