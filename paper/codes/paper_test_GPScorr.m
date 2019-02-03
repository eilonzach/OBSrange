% Do some calculations to test how significant GPS correction is
addpath('/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/OBSrange_v1_MATLAB/functions/');

% datafil = 'SynthBoot_PACMAN_rad1.00_z5000m_fr10_OUT';
folder = 'OUT_OBSrange_synthsurveys_REVISION1_noGPScorr';

dirs = dir(['../figs/figdata/PacificORCA_SynthBoot_surveys_noTAT_REVISION1_table1tests/',folder,'/OUT_wcorr/mats/*.mat']);
for ifil = 1:length(dirs)
    load([dirs(ifil).folder,'/',dirs(ifil).name]);
    fprintf('\n\n::::::::::::\n%s\n',dirs(ifil).name);
    
    % Average shift in x (+ means east of true)
    fprintf('Mean shift in x: %.2f m\n',mean([data.misfit_xsta]));
    
    % Average shift in y (+ means north of true)
    fprintf('Mean shift in y: %.2f m\n',mean([data.misfit_ysta]));
    
    % Average shift in z (+ means station is shallower than true)
    fprintf('Mean shift in z: %.2f m\n',mean([data.misfit_zsta]));
    
    % Average shift in Vw (+ means faster than 1500 m/s)
    fprintf('Mean shift in Vp: %.2f m/s\n',mean([data.misfit_Vw]));
    
    % Calculate average difference in distance of GPS and transponder
    % Account for GPS-transponder offset
    try
        survcog = atan2d(data(1).v_surv_true(:,2),data(1).v_surv_true(:,1));
        [dx,dy] = GPS_transp_correction(data(1).TG_dforward,data(1).TG_dstarboard,survcog');
        [ x_ship_GPS, y_ship_GPS ] = lonlat2xy_nomap( data(1).drop(2), data(1).drop(1), data(1).survlons, data(1).survlats );
        x_ship_TR = x_ship_GPS + dx;
        y_ship_TR = y_ship_GPS + dy;

        r_TR = sqrt( (data(1).obs_loc_xyz(1)*1000-x_ship_TR).^2 + (data(1).obs_loc_xyz(2)*1000-y_ship_TR).^2 );
        r_GPS = sqrt( (data(1).obs_loc_xyz(1)*1000-x_ship_GPS).^2 + (data(1).obs_loc_xyz(2)*1000-y_ship_GPS).^2 );
        dr = mean(r_TR - r_GPS);
        
        % (+ means transponder further away than GPS)
        fprintf('Mean difference in GPS-transponder distance: %.2f m\n',dr);
        
        figure(1); clf;
        subplot(1,2,1)
        plot(x_ship_TR,y_ship_TR,'-ok'); hold on;
        plot(x_ship_GPS,y_ship_GPS,'-or')
%         plot(data(1).survx*1000,data(1).survy*1000,'-b');
        axis equal;
        subplot(1,2,2)
        plot(r_TR - r_GPS); 
    catch
        fprintf('Mean difference in GPS-transponder distance: %.2f m\n',0);
    end
%     pause;
end