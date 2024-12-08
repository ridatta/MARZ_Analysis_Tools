function out =  processBdotSignalsM4(probe_id,block)

 % (1) Load oppositely-polarized signals
 
 load('marz_shot4_bdotprobes.mat'); % CHANGE THIS
 
 % probe_id = 'ABADV05T'; % Probe ID <block>BADV<RadialPosition><Top/Bottom>
 probe_top = [block, probe_id, 'T']; % Top probe
 probe_bottom = [block, probe_id, 'B']; % Bottom Probe
 
 mult = 1;

 
 for ii = 1:length(all_probes)
     % probes{ii} = [all_probes{ii}(end),all_probes{ii}(1:7)]; 
     if strcmpi(polarity{ii},'up')
         p = 'T';
     else
         p = 'B';
     end
     probes{ii} = [all_probes{ii}(5), sprintf('BADV%02d',radial_dist(ii)), p];
 end
 
 top_idx = find(strcmpi(probes,probe_top)); % index of top probe
 bottom_idx = find(strcmpi(probes,probe_bottom)); % index of bottom probe
 
 if strcmpi(string(polarity(top_idx)),'Up')
     probe_top = [probe_top, 'P'];
 else
      probe_top = [probe_top, 'N'];
 end
 
  if strcmpi(string(polarity(bottom_idx)),'Up')
     probe_bottom = [probe_bottom, 'P'];
 else
      probe_bottom = [probe_bottom, 'N'];
 end
 
 t0 = 2800e-9; % start time [s]
 tlm = [0,800e-9]; % integration window, [s]
 k_top = cal_fac(top_idx); % Ts^{-1} / V, calibration factor or 1/A_eff
 k_bottom = cal_fac(bottom_idx); % Ts^{-1} / V
 
 out = getInductiveSignal(time{top_idx},time{bottom_idx},t0,mult*V{top_idx},V{bottom_idx},k_top,k_bottom,tlm);
 out.R = radial_dist(top_idx);
  out.id = [block,probe_id];
 
 showInductiveSignal(out,probe_top,probe_bottom);
 
sgtitle([out.id, ' ' num2str(radial_dist(top_idx)) 'mm'],'FontSize',24);

set(gcf,'Position',[0 0 0.75 1]*1200);
 

end

function out = getInductiveSignal(t1,t2,t0,V1,V2,k1,k2,tlm)

t1 = t1 - t0; % time, [ns]
t2 = t2 - t0; % time, [ns]
V1 = V1 - mean(V1(t1 < 0)); % Zero the voltage
V2 = V2 - mean(V2(t2 < 0)); 

% shift and recenter signals
V1 = V1(t1 >= -500e-9 & t1 <= 2000e-9);
V2 = V2(t2 >= -500e-9 & t2 <= 2000e-9);
if numel(V2) > numel(V1)
    V2 = V2(1:numel(V1)); % ensure size compatibility
    t = t1(t1 >= -500e-9 & t1 <= 2000e-9);
else
    V1 = V1(1:numel(V2));
    t = t2(t2 >= -500e-9 & t2 <= 2000e-9);
end


u1 = -0.5*((V1 - V2)); % signal ~ dPhi/dt

% get stray signal

dtb = (V1 - V2) / (1/k1 + 1/k2); % [T/s]
Vs1 = V1 - dtb / k1; Vs2 = V2 + dtb / k2; % stray voltage [V]

if mean(u1(t > 300e-9 & t < 305e-9)) < 0
    u1 = u1 * -1;
end  

bdot1 = V1 * k1;
bdot2 = V2 * k2; 
bdot_eff = u1 .* 2 ./ (1/k1 + 1/k2); % Inductive signal B_tot [T/s]

% magnetic field
    idx1 = (t >= tlm(1)) & (t <= tlm(end));
    dtB = bdot_eff(idx1);
    out.B = cumtrapz(dtB) * mean(diff(t)); % [T]
    out.tidx = t(idx1); % [s], indices for field
     

out.t = t;
out.V1 = V1 ; 
out.V2 = V2 ; 
out.u1 = u1; % signal ~ dPhi/dt
out.bdot1 = bdot1;
out.bdot2 = bdot2;
out.bdot_eff = bdot_eff;
out.VS = [Vs1, Vs2]; 
end

function showInductiveSignal(out,probe_top,probe_bottom)

t = out.t;
V1 = out.V1; 
V2 = out.V2;
u1 = out.u1;
 
 load('marz_shot3_bdotprobes.mat'); % CHANGE THIS
 
 figure
 subplot(2,2,1)
 plot(t*1e9,V1,'Linewidth',2,'Color',sqclr('b',3),...
        'DisplayName', ['V_1 ', probe_top] ); hold on;
    
        plot(t*1e9,V2,'Linewidth',2,'Color',sqclr('r',5),...
        'DisplayName', ['V_2 ', probe_bottom] ); hold on;
    
     plot(t*1e9,u1,'Linewidth',2,'Color',[0,0,0,0.6],...
        'DisplayName','0.5\times(V_2-V_1)'); hold on;
    

    legend('location','northwest','FontSize',12,'box','off');
    xlabel('time from current start (ns)'); ylabel('Voltage (V)'); 
    xlim([-50,800]); grid on; ylim([-1,1]*700);
    title('Inductive Voltage','FontSize',14); 
    
  subplot(2,2,2)
 plot(t*1e9,out.VS(:,1),'Linewidth',2,'Color',sqclr('b',3),...
        'DisplayName', ['V_{s1} ', probe_top] ); hold on;
    
        plot(t*1e9,out.VS(:,2),'--','Linewidth',2,'Color',sqclr('r',5),...
        'DisplayName', ['V_{s2} ', probe_bottom] ); hold on;
    
    

    legend('location','northwest','FontSize',12,'box','off');
    xlabel('time from current start (ns)'); ylabel('Voltage (V)'); 
    xlim([-50,800]); grid on; ylim([-1,1]*700);
    title('Stray Voltage','FontSize',14);   
 
     subplot(2,2,3)   
  
   plot(t*1e9,out.bdot1/1e9,'Linewidth',2,'Color',sqclr('b',3),...
        'DisplayName', '$\dot{B}_1$'); hold on;
    
   plot(t*1e9,1*out.bdot2/1e9,'Linewidth',2,'Color',sqclr('r',5),...
        'DisplayName', '$\dot{B}_2$'); hold on;
    
   plot(t*1e9,out.bdot_eff/1e9,'Linewidth',2,'Color',[0,0,0,0.6],...
        'DisplayName','$\dot{B}$'); hold on; 
    

    legend('location','northeast','FontSize',14,'Interpreter','latex','box','off');
    xlabel('time from current start (ns)'); ylabel('\partial_t B (T ns^{-1})'); 
    xlim([-50,800]); ylim([-1,1]*0.5); grid on; 
    title('Magnetic Field Rate','FontSize',14); 
    

    subplot(2,2,4) 
     plot(out.tidx*1e9,out.B,'Linewidth',2,'Color','k',...
        'DisplayName', probe_top(1:end-2)); hold on;
       legend('location','northeast','FontSize',14,'Interpreter','latex','box','off');
    xlabel('time from current start (ns)'); ylabel('Magnetic Field  (T)'); 
    xlim([-50,800]); ylim([-1,1]*40); grid on; 
    title('Magnetic Field','FontSize',14); 
end
