
clc; close all; clear;

% This Section loads and shows the raw volatge data and integrated signals
% from all probes. Then saves the data into MAT files for further
% processing.
% The probes should be listed in the EXCEL sheet 'meta_fname'

inDir = checkDir('/Users/rishabhdatta/Dropbox (MIT)/PUFFIN group/Data/Z/MARZ/Shots/MARZ4_z3978/Inductive Probes/');
saveDir = checkDir('/Users/rishabhdatta/Dropbox (MIT)/PUFFIN/Data/MARZ/bdot/z3978/processed/');
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end
meta_fname = checkDir('/Users/rishabhdatta/Dropbox (MIT)/PUFFIN group/Data/Z/MARZ/Shots/MARZ4_z3978/Inductive Probes/InductiveProbeLabelsMARZ copy.xlsx'); % Metadata

meta_data = importfileMARZ4(meta_fname); % meta data of probes

all_probes = meta_data.SignalName; % all probe names
cal_fac = meta_data.CalFactor; % all probe names
probe_num = meta_data.Label; % probe numbers
radial_dist = meta_data.Radius; % [mm]
polarity = meta_data.Polarity; % polarity


% sort by distance
[radial_dist,idx] = sort(radial_dist);
all_probes = all_probes(idx);
cal_fac = cal_fac(idx);
probe_num = probe_num(idx);
polarity = cellstr(polarity); polarity = polarity(idx);

tidx = 2800+[0,800];

all_data = importMARZ4Data(checkDir('/Users/rishabhdatta/Dropbox (MIT)/PUFFIN group/Data/Z/MARZ/Shots/MARZ4_z3978/Inductive Probes/InductiveProbes.txt')); % Raw data

dist = 10;

pidx = find(radial_dist <= dist); % only plot data <= dist
% pidx = 1:length(all_probes);

f1 = figure(); f2 = figure(); f3 = figure();
% Iterate over each probe
for kk = 1:numel(pidx)
    inum = pidx(kk);
    probe_name = upper(all_probes{inum});
    
    str_time = ['Timens' probe_name];
    str_volt = ['Amplitude' probe_name];

    t = all_data{:,str_time} * 1e-9; time{inum} = t; % ns 
    V{inum} = all_data{:,str_volt}; % V
    dtB{inum} =   V{inum} .* cal_fac(inum); % [T/s]
    dtB{inum} = dtB{inum} - mean(dtB{inum}(t < 2800e-9));
    

    % integrate

    idx1 = (t >= tidx(1)*1e-9) & (t <= tidx(end)*1e-9);
    bdot = dtB{inum};
    bdot = bdot(idx1);
    B{inum} = cumtrapz(bdot) * max(diff(t)); % [T]
    
   lsty = {'-',':','-.','--','-',':','-.','--',};
   clr = lines(6);
   
   if radial_dist(inum) <= 30
   colr = clr( (radial_dist(inum)-5) / 5 + 1, :); 
   else
       colr = clr( (radial_dist(inum)-20) / 5 + 1, :); 
   end
% colr = clr(kk,:);

    figure(f1);
    plot(t*1e9,dtB{inum},'LineStyle',lsty{kk},'Linewidth',1,'Color',colr,...
        'DisplayName',[num2str(radial_dist(inum)) 'mm (', ...
        char(polarity(inum)), ') ' upper(strrep(probe_name,'_','\_'))]); hold on;
    
    figure(f2);
    plot(t*1e9,V{inum},'LineStyle',lsty{kk},'Linewidth',1,'Color',colr,...
        'DisplayName',[num2str(radial_dist(inum)) 'mm (', ...
        char(polarity(inum)), ') ' upper(strrep(probe_name,'_','\_'))]); hold on;
    
      figure(f3);
    plot(t(idx1)*1e9,B{inum},'LineStyle',lsty{kk},'Linewidth',1,'Color',colr,...
        'DisplayName',[num2str(radial_dist(inum)) 'mm (', ...
        char(polarity(inum)), ') ' upper(strrep(probe_name,'_','\_'))]); hold on;
    
  
    end

figure(f1);
formatPlots(1200); 
    legend('location','northwest','FontSize',14,'NumColumns',2);
    set(gcf,'Position',[0 0 2 1]*500);
    xlabel('t (ns)'); ylabel('\partial_t B (T/s)'); 
    xlim([0,500]+2800); grid on; figure(f1);
    title('\partial_t B');
    saveas(gcf,[saveDir,'bdot-dist-' num2str(dist) '.png']); 
    
    
 figure(f2);   
    formatPlots(1200); 
    set(gcf,'Position',[0 0 2 1]*500);
    legend('location','northwest','FontSize',14,'NumColumns',2);
    xlabel('t (ns)'); ylabel('Voltage (V)'); 
    xlim([0,500]+2800); grid on; ylim([-1,1]*700);
       title('Raw Voltage');
        xl = xticklabels();
    xticklabels(str2num(cell2mat(xl))-2800);
    saveas(gcf,[saveDir,'V-dist-' num2str(dist) '.png']); 

    figure(f3);  
    formatPlots(1200); 
   set(gcf,'Position',[0 0 2 1]*500);
    legend('location','northwest','FontSize',14,'NumColumns',2);
    xlabel('t (ns)'); ylabel('B (T)'); 
    xlim([0,500]+2800); grid on; ylim([-1,1]*35);
    xl = xticklabels();
    xticklabels(str2num(cell2mat(xl))-2800);
    
    title('B (T)');
    saveas(gcf,[saveDir,'B-dist-' num2str(dist) '.png']); 


    % Rename the misnamed probes

   idx =  find(strcmpi(all_probes,"BADCV7"));
   all_probes{idx} = 'BADVC7';

   idx = find(strcmpi(all_probes,"BADCV8"));
   all_probes{idx} = 'BADVC8';

    
    save('marz_shot4_bdotprobes.mat','meta_data','all_probes',...
        'cal_fac',...
        'probe_num',...
        'polarity',...
        'radial_dist',...
        'time','V','dtB','B'); % save as MAT file
% %     
%%


% This section combines opposite polarity signals to determine the
% inductive and stray signals. Then saves the processed data for further
% analysis. It also plots all processed signals for V_ind, B_dot, and B.

% NOTE: This section is incomplete. processBdotSignalsM4.m must be modified
% to be compatible with the new naming scheme used in MARZ 4.

clc; close all; clear;
 
load('marz_shot4_bdotprobes.mat');

saveDir = checkDir('/Users/rishabhdatta/Dropbox (MIT)/PUFFIN/Data/MARZ/bdot/z3978/processed/');
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end

xlm = [0,400];
d = [5, 10];

% Block C
for ii =1:length(d)
    
        probe_id = sprintf('BADV%02d',d(ii));
        out = processBdotSignalsM4(probe_id,'C'); % process signal using custom function for MARZ 4 signals
        save([saveDir,'r=' num2str(out.R) 'mm_' out.id  '.mat'],'out');
        saveas(gcf, [saveDir, 'r=' num2str(out.R) 'mm_' out.id '.png']); 

    close all;
end

% Block A - note on block A 2 probes failed
for ii =1:length(d)
    
        probe_id = sprintf('BADV%02d',d(ii));
        out = processBdotSignalsM4(probe_id,'A'); % process signal using custom function for MARZ 4 signals
        save([saveDir,'r=' num2str(out.R) 'mm_' out.id  '.mat'],'out');
        saveas(gcf, [saveDir, 'r=' num2str(out.R) 'mm_' out.id '.png']); 

    close all;
end



% Plot all data on same figure


figure % Voltage
% for ii =1:length(d)
%         probe_id = sprintf('ABADV%02d',d(ii));
%         load([saveDir,'r=' num2str(d(ii)) 'mm_' probe_id '.mat'],'out');
%         plot(out.t*1e9,(out.u1),'LineWidth',2.5,'Color',sqclr('r',ii),...
%             'DisplayName',[num2str(out.R) 'mm ' probe_id]); hold on;
% end

for ii =1:length(d)
        probe_id = sprintf('CBADV%02d',d(ii));
        load([saveDir,'r=' num2str(d(ii)) 'mm_' probe_id '.mat'],'out');
        plot(out.t*1e9,(out.u1),'LineWidth',2.5,'Color',sqclr('b',ii),...
            'DisplayName',[num2str(out.R) 'mm ' probe_id]); hold on;
end


 set(gcf,'Position',[0 0 1 0.6]*1000);
 formatPlots(700,2.5); legend('location','northwest','NumColumns',2); 
    xlabel('time from current start (ns)'); ylabel('Voltage (V)'); 
    xlim(xlm); grid on; ylim([-300,600]);
    title('Inductive Voltage'); legend('location','northwest');
    saveas(gcf, [saveDir, 'bdot_voltage.png']); 
    
f1 = figure; % B-field
% for ii =1:length(d)
%         probe_id  = sprintf('ABADV%02d',d(ii));
%         load([saveDir,'r=' num2str(d(ii)) 'mm_' probe_id '.mat'],'out');
%         plot(out.tidx*1e9,(out.B),'LineWidth',2.5,'Color',sqclr('r',ii),...
%             'DisplayName',[num2str(out.R) 'mm ' probe_id]); hold on;
%         % save as csv
%         writematrix([out.tidx,out.B],[saveDir, 'r=' num2str(d(ii)) 'mm_' probe_id '.csv']) 
% end
for ii =1:length(d)
        probe_id  = sprintf('CBADV%02d',d(ii));
        load([saveDir,'r=' num2str(d(ii)) 'mm_' probe_id '.mat'],'out');
        plot(out.tidx*1e9,(out.B),'LineWidth',2.5,'Color',sqclr('b',ii),...
            'DisplayName',[num2str(out.R) 'mm ' probe_id]); hold on;
        % save as csv
        writematrix([out.tidx,out.B],[saveDir, 'r=' num2str(d(ii)) 'mm_' probe_id '.csv']) 
end    
 

% T probe
%  plot(tp.t(bidx)*1e9,tp.B,'LineWidth',2.5,'Color',sqclr('r',ii+1),...
%             'DisplayName',[num2str(tp.R) 'mm ' tp.id]); hold on;
        
 %        tx = linspace(0,600,100);
 % T = 600; I = 20 .* (sin(1*pi*tx/T)).^2;       

 set(gcf,'Position',[0 0 1 0.6]*1000);
 formatPlots(700,2.5); legend('location','northwest','NumColumns',2); 
    xlabel('time from current start (ns)'); ylabel('B-field (T)'); 
    xlim(xlm); grid on; ylim([0,50]);
    title('Magnetic Field'); legend('location','northwest'); 
    
    
%     yyaxis right
%     plot(tx,I,'LineWidth',3,'Color',[0.5 0.5 0.5 0.75],'DisplayName','Current'); 
%     ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.YAxis(2).Color = 'k';
% 
%        ylim([0,30]); ylabel('Current (MA)');
     saveas(gcf, [saveDir, 'bfield.png']); 
     saveas(gcf, [saveDir, 'bfield-aps-2.fig']); 
     
     
     
     
 figure % BStray-Voltage
% for ii =1:length(d)
%     
%         probe_id  = sprintf('ABADV%02d',d(ii));
%         load([saveDir,'r=' num2str(d(ii)) 'mm_' probe_id '.mat'],'out');
%         plot(out.t*1e9,(out.VS(:,2)),'LineWidth',2.5,'Color',sqclr('r',ii),...
%             'DisplayName',[num2str(out.R) 'mm ' probe_id]); hold on;
%     
% end

for ii =1:length(d)
    
        probe_id  = sprintf('CBADV%02d',d(ii));
        load([saveDir,'r=' num2str(d(ii)) 'mm_' probe_id '.mat'],'out');
        plot(out.t*1e9,(out.VS(:,2)),'LineWidth',2.5,'Color',sqclr('b',ii),...
            'DisplayName',[num2str(out.R) 'mm ' probe_id]); hold on;
    
end

 set(gcf,'Position',[0 0 1 0.6]*1000);
 formatPlots(700,2.5); legend('location','northwest','NumColumns',2); 
    xlabel('time from current start (ns)'); ylabel('Voltage (V)'); 
    xlim(xlm); grid on; ylim([-300,600]*1);
    title('Electrostatic Voltage'); legend('location','northwest');  
     saveas(gcf, [saveDir, 'stray_voltage-2.png']);   
     
     
%% Velocity

% This section caluclates the velocity from time of flight of B-dot
% signals. Note the previous two sections must be run before running this.
% Point selecttion is manual. Set new = true if you want to re-do points
% selection.

clc; close all; clear;

saveDir = checkDir('/Users/rishabhdatta/Dropbox (MIT)/PUFFIN/Data/MARZ/bdot/z3978/processed/'); % save path
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end

% Specify the two probe signals to compare
 tag = '-5mm-10mm-1'; % file tag, description

 load([saveDir,'r=5mm_CBADV05.mat'],'out'); % Load first probe, see .mat files in saveDir
 loc1 = out;
 
 load([saveDir,'r=10mm_CBADV10.mat'],'out'); % Load second probe, see .mat files in saveDir
 loc2 = out;
 xlm= [0,500];
  

% find peaks in signal

figure % plot the probe signals
set(0, 'DefaultFigureRenderer', 'painters');
h1 = plot(loc1.t*1e9,loc1.u1,'LineWidth',2,'Color',sqclr('r',3),...
    'DisplayName',[num2str(loc1.R) 'mm ' loc1.id]); hold on;
h2 = plot(loc2.t*1e9,loc2.u1,'LineWidth',2,'Color',sqclr('r',6),...
    'DisplayName',[num2str(loc2.R) 'mm ' loc2.id]); hold on;

xlabel('time (ns)'); ylabel('Signal (V)'); 
 set(gcf,'Position',[0 0 1 0.75]*1800);
    xlabel('time form current start (ns)'); ylabel('Voltage (V)'); 
 grid on; ylim([-200,700]);
     xlim(xlm); 
    formatPlots(900,2.5);
    title('Voltage'); legend('location','northwest');


new = false; % set to TRUE if you want to re-do points selection

% Point selectionis done manually
% first, select points on first probe then double click
% next, selct points on sencd signal then doible click
if new 
h1.Color(4) = 1; h2.Color(4) = 0.3;
[x1,y1] = getpts; % this code gets the newpoints
x1 = x1(1:end-1);
y1 = y1(1:end-1);
plot(x1,y1,'O','Color',sqclr('r',3),'MarkerSize',25,'HandleVisibility','off','Linewidth',2); 

h1.Color(4) = 0.5; h2.Color(4) = 1;
[x2,y2] = getpts;
x2 = x2(1:end-1); 
 y2 = y2(1:end-1); 
plot(x2,y2,'O','Color',sqclr('r',6),'MarkerSize',25,'HandleVisibility','off','Linewidth',2); 

save([saveDir, tag 'velocity_analysis.mat'],'x1','y1','x2','y2');

else
    load([saveDir, tag 'velocity_analysis.mat'],'x1','y1','x2','y2');
    h1.Color(4) = 1; h2.Color(4) = 1;
    plot(x1,y1,'O','Color',sqclr('r',3),'MarkerSize',30,'HandleVisibility','off','Linewidth',3); 
    plot(x2,y2,'O','Color',sqclr('r',6),'MarkerSize',30,'HandleVisibility','off','Linewidth',3); 
end

 h1.Color(4) = 1; h2.Color(4) = 1;


% connect the selected points
for ii = 1:length(x1)
    plot([x1(ii),x2(ii)],[y1(ii),y2(ii)],'Color','k','LineWidth',2,'HandleVisibility','off'); hold on;
end

title(tag(2:end)); 
set(gcf,'Position',[0 0 1000 500]);
saveas(gcf,[saveDir tag 'trasit_time.png']);

       % SAVE FOR POP
    formatPlots(900,2.5); 
    title('');
    % yticks([0,300,600]);
    ylabel('Voltage [V]')
    % set(gca,'FontSize',44)
    xlim([0,400]);
    ylim([0,300]);
    legend('off');
    savePath = checkDir('/Users/rishabhdatta/Dropbox (MIT)/PUFFIN/Data/MARZ/Figures_POP/');
    saveas(gcf,[savePath,'z3781_bdot_transitTime.png']);


% check how good the point selection is by caluclatibg derivative of the
% signals

figure
sm_u1 = smoothdata(loc1.u1,'movmedian',20); sm_u2 = smoothdata(loc2.u1,'movmedian',20);
subplot(2,1,1);
plot(loc1.t*1e9,loc1.u1,'Color',[sqclr('r',3),0.7]); hold on; % first probe
plot(loc1.t*1e9,sm_u1,'Color','k'); hold on;

plot(loc2.t*1e9,loc2.u1,'Color',[sqclr('r',6),0.6]); hold on; % second probe
plot(loc2.t*1e9,sm_u2,'Color','k'); hold on;

xlim([0,800]); xlabel('t (ns)'); ylabel('Voltage (V)'); 

for ii = 1:numel(x2)
plot([x1(ii),x1(ii)],[-1,1]*700,'Color',sqclr('r',3),'HandleVisibility','off'); hold on;
plot([x2(ii),x2(ii)],[-1,1]*700,'r','HandleVisibility','off'); hold on;
end
legend('off');
ylim([-700,700]);


plot(x1,y1,'bo','MarkerSize',15); formatPlots(); grid on;
plot(x2,y2,'ko','MarkerSize',15); formatPlots(); grid on;
formatPlots(); legend('off');

% Now plot time-deribative of smoothed signal

subplot(2,1,2);
plot(loc1.t(2:end)*1e9,diff(sm_u1)./diff(loc1.t),'Color',sqclr('r',3));  hold on;
plot(loc2.t(2:end)*1e9,diff(sm_u2)./diff(loc2.t),'Color',sqclr('r',6));  hold on;

for ii = 1:numel(x2)
plot([x1(ii),x1(ii)],[-1,1]*4e10,'Color',sqclr('r',3)); hold on;
plot([x2(ii),x2(ii)],[-1,1]*4e10,'r'); hold on;
end
formatPlots(); grid on;
xlim([0,800]); 
ylim([-1,1]*4e10); legend('off');
set(gcf,'Position',[0 0 1 1]*900);

% VELOCITY CALCULATION

sort(x1); sort(x2); 
delt = x2 - x1; % transit time [ns]
R = 5e-3; % [m] SPECIFY PROBE SEPARATION
V = R ./ (delt * 1e-9); % m/s

% Uncertainty
delR = 1e-3; % uncertainty in space
deldelt = sqrt(2^2 + 2^2) * 1e-9; % uncertainty in time
% deldelt = 0.1*delt*1e-9;
delV = sqrt( (-R./(delt*1e-9).^2).^2 .* deldelt.^2 + ...
    (1./ (delt*1e-9)).^2 .* delR.^2 * 1); % voltage uncertianty; error propgation [m/s]


xm = 0.5*(x1+x2); % time at midpoint, ns
% Plot at midpoint
figure
plot(xm,V*1e-3,'O','Color','k','Linewidth',4,'MarkerSize',15,'HandleVisibility','off'); hold on;
errorbar(xm,V*1e-3,delV*1e-3,delV*1e-3,(x2-x1)/2,(x2-x1)/2,'O',...
    'Linewidth',2,'MarkerSize',15,...
    'Color','k','MarkerEdgeColor','k',...
    'DisplayName','Average Velocity'); hold on;

% Plot Linear fit

[PP,S] = polyfit(xm,V,1);
xx = linspace(min(xm),max(xm),100)';
[V_fit,del] = polyval(PP,xx,S);

xx2 = linspace(min(xm),500,100);
[V_fit2,del2] = polyval(PP,xx2,S);

s = shadedErrorBar(xx,V_fit*1e-3,0*del*1e-3,...
    'lineProps',{'-','Color',[sqclr('b',2),1],'Linewidth',4,'DisplayName','Fit'}); 
s = shadedErrorBar(xx2,V_fit2*1e-3,del*1e-3,...
    'lineProps',{'--','Color',[sqclr('b',2),1],'Linewidth',4,'HandleVisibility','off'}); 

formatPlots();



save([saveDir, tag 'velocity_analysis.mat'],'x1','y1','x2','y2','V','delV');

xlabel('time from current start (ns)'); ylabel('Velocity (km/s)');
 ylim([0,400]); 
 grid on; set(gcf,'Position',[0 0 1200 1200/1.8]);
title(tag(2:end)); 
    xlim(xlm); 
    formatPlots(900,2.5);

set(gcf,'Position',[0 0 1000 500]);
saveas(gcf,[saveDir tag 'velocity_displacement.png']);

         % SAVE FOR POP
    formatPlots(900,2.5); 
    title('');
    xlabel('t [ns]')
    ylabel('V [km/s]')
    % set(gca,'FontSize',44)
    % xlim([0,500]);
    xlim([0,400]);
    ylim([0,300]);
    legend('off');
%     savePath = checkDir('/Users/rishabhdatta/Dropbox (MIT)/PUFFIN/Data/MARZ/Figures_POP/');
%     saveas(gcf,[savePath,'z3781_bdot_velocity.png']);
% 
% plot(data(:,1),data(:,2)/1e3,'--r','Linewidth',3,Displayname='Simulation');


     