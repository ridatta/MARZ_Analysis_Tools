classdef SVSAnalysis
    properties
    end
    methods
        function obj = SVSAnalysis()            
        end
        function out = loadRawImage(~,svsobj)
            % out = image matrix
            % svsobj = svsobj with fields inDir [inpput dir] and
            % fid [filename]
            out =  double(hdfread([svsobj.inDir, svsobj.fid, '.hdf'],'fore'));
        end
        function out = correctDistortion(~,svsobj)
            % Corrects for strek image distrotion
            % Requires the SMASH toolbox from Sandia National Labs
            import SMASH.ImageAnalysis.Image
            import SMASH.SpectroscopyAnalysis.SVSAnalysis.Spectroscopy

            data=Image([svsobj.inDir, svsobj.fid, '.hdf'],'sydor');
            obj=Spectroscopy(data,svsobj.inum,'streak','Auto','False');
            % apply distortion correction for the 500 ns sweepts
            obj=obj.distortion('speed','500ns');
            obj.distortionCal=obj.UnmodifiedData.Data;

            out = obj.ActiveImage.Data;
        end
        function showSVSImage(~,svsobj,field,varargin)
            % Shows the image stored in a given field
            figure
            imagesc(svsobj.(field)/1e3);
            set(gca,'Ydir','Normal');
            colormap(hot)
            title( [strrep(svsobj.fid,'_','\_'), ' ',  strrep(field,'_','\_')] );
            if isempty(varargin)
                ylabel('wavelength (a.u.)');
                xlabel('time (a.u.)');
            else
                xlim(svsobj.calib.t2px([2850,3300]))
                xlabels = 2800:100:3400; % ns
                ylabels = 400:100:700; % nm
                  xticks(svsobj.calib.t2px(xlabels'));
%                   xticklabels(xlabels'-2800);
                  xticklabels(xlabels'-0);
                  yticks(svsobj.calib.wl2px(ylabels'));
                  yticklabels(ylabels');
                   
                  ylim(svsobj.calib.wl2px([400,700]));
                 ylabel('wavelength (nm)');
                 xlabel('time (ns)'); 
            end
            set(gcf,'Position', [262   118   820   556]);
            set(gca,'TickDir','out')
            cb = colorbar(); ylabel(cb,'counts (\times 10^3 a.u.)')
            formatPlots();
            
            % save
            exportgraphics(gcf,[svsobj.saveDir, svsobj.fid '_raw.tiff']);
        end
        function [laser_img,laser_img_corr] = getLaserCalibImg(obj,fid,svsobj)
            % loads the laser clibration prehot images
            % fid refers to the laser image filename
            l = cell(1,numel(fid));
             for ii = 1:length(fid)
                l{ii} = struct('inum',5,'fid',fid(ii),'inDir',svsobj.inDir,'saveDir',svsobj.saveDir);
                l{ii}.img = obj.loadRawImage(l{ii});
                l{ii}.img_corr = obj.correctDistortion(l{ii});
            end

        % create combined image
        laser_img = zeros(size(l{1}.img));
        laser_img_corr = zeros(size(l{1}.img));
        for ii = 1:length(fid)
                laser_img = laser_img + l{ii}.img;
                laser_img_corr = laser_img_corr + l{ii}.img_corr;
        end
        end
        function checkAlignment(~,svsobj,field,tid)
            % Checks alignment by taking lineouts at times tid [px] of image
            % stored in 'filed' of svsobj 
            [m,~] = size(svsobj.img);
            for ii = 1:length(tid)
                out{ii} = sum(svsobj.(field)(:,tid(ii)-10:tid(ii)+10),2)/21;
            end

            figure
            for ii = 1:length(tid)
                plot(1:m,out{ii},'Color',sqclr('b',ii),'DisplayName',['t = ' num2str(tid(ii)) ' a.u. +/- 20 a.u.'],'LineWidth',1.5); hold on;
            end
            xlabel('Wavelength (a.u.)'); ylabel('Counts (a.u.)'); 
            xlim([700,1050]); grid on;
            formatPlots(600); set(gcf,'Position',[0 0 1200 600]); legend('NumColumns',2,'Location',"best")
        end
        function out = wlcalib(obj,svsobj,field,wl,tidx)
            % Perform wavelength calibration by comparing position of laser lines
            % with their known wavlengths.
            % wl = known wavelengths of laser lines [nm]
            % tidx = [px] time index to take lineouts at
             [m,~] = size(svsobj.img); Lq = 100:m-100;
             Iq = sum(svsobj.(field)(Lq,tidx-20:tidx+20),2)/41; % get intenisty lineout
             Iq = Iq - mean(Iq(1:10)); % zero signal
             % find location of peaks i.e. laser lines
             % depending on the laser line strengthm you may have to
             % change the 'MinPeakHeight'
             [~,idx] = findpeaks(Iq,'MinPeakHeight',max(Iq)/5,'NPeaks',numel(wl)); 
             figure();
             for ii = 1:length(wl) % fit gaussians to the laser lines
                 [sig(ii),mu(ii),A(ii)] = mygaussfit(Lq(idx(ii)-20:idx(ii)+20),Iq(idx(ii)-20:idx(ii)+20));
                 subplot(1,length(wl),ii)
                 plot(Lq(idx(ii)-20:idx(ii)+20),Iq(idx(ii)-20:idx(ii)+20),'o','Markersize',10); hold on;
                 x = linspace(min(Lq),max(Lq),20*numel(Lq)); 
                 plot(x,obj.gaussianFn(x,mu(ii),sig(ii),A(ii)),'-'); hold on;
                 xlabel('Wavelength (a.u.)'); ylabel('Counts (a.u.)'); formatPlots();
                 xlim([Lq(idx(ii)-20),Lq(idx(ii)+20)]);
             end
             set(gcf,'Position',[0 0 1200 600]);
             
             
             P = polyfit(mu,wl,1); % Px to wl [nm]
             
            figure
            plot(mu,wl,'o',"markerSize",20,'LineWidth',2,'HandleVisibility','off'); hold on;
            str = sprintf('WL [nm] = %0.3f * px + %0.1f',P(1),P(2));
            xx = linspace(0,m,100); yy = polyval(P,xx);
            plot(xx,yy,'--','LineWidth',2,'DisplayName','fit');
            text(mean(xx)+100,mean(yy), str,'FontSize',24);
            ylabel('Wavelength (nm)'); xlabel('y-position (px)')
            grid on; formatPlots();
            xlim([0,m])
            set(gcf,'Position', [0   0   1200   600]*1.5);
            
            % return calibration functions
            out.px2wl = @(x) P(1) * x + P(2) 
            out.wl2px = @(wl) (wl - P(2)) / P(1)
        end
        function out = gaussianFn(~,x,mu,sigma,A) % returns a gaussian
            out = A .* exp( -1 * (x - mu).^2 / (2 * sigma.^2)); 
        end
        function out = timeCalib(~,svsobj,field,delt,machine_time,window)
            % Performs time calibration based on streak image stored in
            % 'field' of svsobj
            % delt [ns] = comb separation
            % machine_time [ns] = machine time
            data = svsobj.(field); % get svs image

            % Get comb
            comb = data(window(1):window(end),:); comb = sum(comb,1); % temporary variable
            [m,n] = size(comb);
            [pks,idx] = findpeaks(comb,'MinPeakDistance',50,'MinPeakHeight',2e5);

            % Show machine time marker

            impulse = data(1700:1900,:);
            impulse = sum(impulse,1);
            [m,n] = size(impulse);
            [mkr,mkr_idx] = findpeaks(impulse,'MinPeakDistance',50,'MinPeakHeight',max(impulse)/4);
            [~,maxidx] = max(mkr);
            mkr = mkr(maxidx); mkr_idx = mkr_idx(maxidx);

            % get dispersion relation gradient
        
            t1 = 0; tend = t1 + (numel(pks)-1) * delt;
            t = linspace(t1,tend,numel(pks));
            P = polyfit(idx,t,1); % px to time [ns]
            m = P(1);

            % correct for machine time
            c = machine_time - m * mkr_idx; 

            out.px2t = @(x) m * x + c; % px to ns
            out.t2px = @(t) (t - c) / m; % ns to px

            % plot
            figure
            subplot(2,2,1);
            plot(1:n,comb,'linewidth',1.5); hold on;
            plot(idx,pks,'om','markersize',15); 
            xlim([0,n]); grid on;
            ylabel('counts (a.u.)');
            formatPlots(); title('Comb');
            xticks(0:500:2000);xticklabels([]);
            legend('off');

            subplot(2,2,2);
            plot(1:n,impulse,'linewidth',1.5); hold on;
            plot(mkr_idx,mkr,'om','markersize',15); 
            xlim([0,n]); grid on;
            ylabel('counts (a.u.)');
            formatPlots(); title('Machine time marker');
            xticks(0:500:2000);xticklabels([]);
            legend('off');

            subplot(2,2,3);
            plot(idx,t,'o','markersize',15); hold on;
            xx = linspace(0,n,100); % px
            yy = polyval(P,xx); % ns
            plot(xx,yy,'-','color','r','linewidth',1.5); 
            xlabel('x (px)');
            xticks(0:500:2000);
            ylabel({'time from', 'first marker (ns)'});
            legend('off');
            str = sprintf('t = %0.2f * px + %0.2f',P(1),P(2));
            text(mean(xx)-500,mean(yy),str,'FontSize',12);
            grid on; formatPlots();
            xlim([0,n])
            legend('off');


            subplot(2,2,4);
            px = 1:n; t = out.px2t(px);
            plot(px,t,'-','color','b','linewidth',1.5); hold on;
            plot(mkr_idx,machine_time,'ro','Markersize',15);
            xlabel('x (px)');
            ylabel('t (ns)');
            legend('off');
            str = sprintf('t = %0.2f * px + %0.2f',m,c)
            text(mean(px)-500,mean(t),str,'FontSize',12);
            grid on; formatPlots();
            xlim([0,n]); xticks(0:500:2000);
            set(gcf,'Position', [0   0   1600   1400]*1.5);
            legend('off');

        end
        
        function [out, zeroval] = zeroCorrect(obj,svsobj,field,rows,cols)
            % Correct for non-zero background
            % rows, cols set window of dark pixels
            img = svsobj.(field)(rows(1):rows(2),cols(1):cols(2));
            zeroval = mean(img(:));
            out = svsobj.(field) - zeroval; % returns zer (dark counts) corrected image
        end
        function fwhm = getInstrumBroadening(obj,svsobj,field,tid)
        % Get instrument broadening base don laser preshot image
        % tid = times to pick [ns]
        dt = 20; % [ns]
        [out,wl] = obj.getLineout(svsobj,field,tid,dt);
        out = out{1};
        % get lineouts
        figure
            out = out - mean(out(1:10)); % zero the signal
            idx{1} = (wl > (458-10)) & (wl < (458+10));
            idx{2} = (wl > (543-10)) & (wl < (543+10));
            spec{1} = out(idx{1}); spec{2} = out(idx{2});
            figure
            for ii = 1:numel(spec)
                subplot(1,2,ii)
                plot(wl(idx{ii}),spec{ii}/max(spec{ii}),'o','MarkerSize',10,'DisplayName',['t = ' num2str(tid) ' ns \pm' num2str(dt)  ' ns']); hold on;
                % get gaussian fit
                [sigma(ii),mu,A]= mygaussfit(wl(idx{ii}),spec{ii}/max(spec{ii}));
                xx = linspace(min(wl(idx{ii})),max(wl(idx{ii})),20 * numel(wl(idx{ii})));
                plot(xx,obj.gaussianFn(xx,mu,sigma(ii),A),'Linewidth',2,'Color','#D95319','DisplayName','Gaussian Fit');
                text(mu,0.5,['FWHM = ' num2str(2.355 * sigma(ii),4) ' nm'],'Fontsize',12);
                xlabel('Wavelength (nm)'); ylabel('Counts (ns)'); 
                grid on; ylim([0,1.4]);
                formatPlots(600); legend('Location',"best");
            end
        set(gcf,'Position',[0 0 1200 600]); 
        fwhm = 2 * sqrt(2 * log(2)) .* sigma;
        end
       
        function [out,wl] = getLineout(~,svsobj,field,tid,dt)
            % Returns lineout at time tid averaged over dt for image
            % stored in 'field' of svsobj
            % tid [ns]
            % dt [+/- ns]
            % returns avg. lineout at tid +/- dt vs. wavelength (px)
            px_tid = svsobj.calib.t2px(tid); px_tid = round(px_tid);
            px_dt = round(svsobj.calib.t2px(1000+dt) - svsobj.calib.t2px(1000));
            for ii = 1:length(tid)
                out{ii} = sum(svsobj.(field)(:,px_tid(ii)-px_dt:px_tid(ii)+px_dt),2);
                out{ii}  = out{ii} / (px_dt * 2 + 1); % avraging
                out{ii} = out{ii} - mean(out{ii}(1:10)); % zero teh signal
            end
            [m,~] = size(svsobj.(field));
            wl = svsobj.calib.px2wl(1:m);
        end
        function [out,wl] = getLineoutCorr(~,svsobj,field,tid,dt)
            % Returns lineout at time tid averaged over dt for image
            % stored in 'field' of svsobj after correction for instrument
            % response
            % tid [ns]
            % dt [+/- ns]
            % returns avg. lineout at tid +/- dt vs. wavelength (px)
            px_tid = svsobj.calib.t2px(tid); px_tid = round(px_tid);
            px_dt = round(svsobj.calib.t2px(1000+dt) - svsobj.calib.t2px(1000));

            [m,~] = size(svsobj.(field)); % wavelength
            wl = svsobj.calib.px2wl(1:m);
            idx = wl >= 400 & wl <=700;
            wl = wl(idx);

            for ii = 1:length(tid)
                out{ii} = sum(svsobj.(field)(:,px_tid(ii)-px_dt:px_tid(ii)+px_dt),2);
                out{ii}  = out{ii} / (px_dt * 2 + 1); % avraging
                out{ii} = out{ii} - mean(out{ii}(1:10)); % zero the signal
                out{ii} = out{ii}(idx); % correct range
                out{ii} = out{ii} .* svsobj.instrum_resp_fn(wl');
            end
            
        end
        function out = t_shift(~,img,shift)
            % img = [mxn] size image
            % shift = [m x 1] size shifts in px
            % shift an image by px values given in shift
            [m,~] = size(img);
            out = 0 * img;
            for ii = 1:m
                a = img(ii,:); 
                b = circshift(a,shift(ii)); 
                if shift(ii) > 0
                    b(1:shift(ii)) = 0;
                else
                    b(end+shift(ii)+1:end) = 0;
                end
                out(ii,:) = b;
            end
        end
        function [a,b,corr] = find_shift_scale(~,x,y_target,y)

            % x = wavelength 
            % y_target = target spectrum or fast spectrum
            % y = slow spectrum


            PP = spline(x,y);
            func = @(x) ppval(PP,x); % fitting func.
            ft = fittype(@(a,b,x) func((x+b)/a),...
                'coefficients',{'a','b'},'independent',{'x'});


            f = fit(x(:),y_target(:),ft,'StartPoint',[1,25]);

            a = f.a; b = f.b; 

            corr = f(x);

            end


            function out = sweep_correct(~,x,y,shift,scale)

            % x = wavelength 
            % y = slow spectrumS
            % shift 
            % scale

                PP = spline(x,y);
                out = ppval(PP,(x+shift)/scale);
            end
    end
end