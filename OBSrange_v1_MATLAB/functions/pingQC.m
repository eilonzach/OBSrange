function [ datao, datao_bad ] = pingQC( data, vp_w, res_thresh )
%
% Ping quality control
%
lat_drop = data.lat_drop;
lon_drop = data.lon_drop;
z_drop = data.z_drop;
lats_ship = data.lats;
lons_ship = data.lons;
twt = data.twt;

[ x_ship, y_ship ] = lonlat2xy( lon_drop, lat_drop, lons_ship, lats_ship );
% [ x_drop, y_drop ] = lonlat2xy( lon_drop, lat_drop, lon_drop, lat_drop );
x_drop=0;
y_drop=0;

twt_pre = calcTWT(x_drop, y_drop, z_drop, 0, 0, x_ship, y_ship, 0, vp_w);
dtwt = twt - twt_pre;

I_bad = abs(dtwt)*1000 > res_thresh;
if (sum(I_bad)/length(I_bad))>0.9
    % MUST BE 1-way travel time??!!
    twt = 2*twt;
    dtwt = twt - twt_pre;
    I_bad = abs(dtwt)*1000 > res_thresh;
end
if (sum(I_bad)/length(I_bad))>0.9
    error('Something very wrong with travel times...');
end

% Plot
if 0
    obs = 1:length(dtwt);
    figure(50); clf;
    plot(obs,dtwt*1000,'ok'); hold on;
%     plot( [obs(1),obs(end)] , [(mean(dtwt)+2*std(dtwt))*1000, (mean(dtwt)+2*std(dtwt))*1000] ,'-');
%     plot( [obs(1),obs(end)] , [(mean(dtwt)-2*std(dtwt))*1000, (mean(dtwt)-2*std(dtwt))*1000] ,'-');
    xlabel('Observation #');
    ylabel('Residual (ms)');
end

datao.lat_drop = lat_drop;
datao.lon_drop = lon_drop;
datao.z_drop = z_drop;
datao.lats = lats_ship(~I_bad);
datao.lons = lons_ship(~I_bad);
datao.t_ship = data.t_ship(~I_bad);
datao.twt = twt(~I_bad);
datao.sta = data.sta;

datao_bad.lats = lats_ship(I_bad);
datao_bad.lons = lons_ship(I_bad);
datao_bad.t_ship = data.t_ship(I_bad);
datao_bad.twt = twt(I_bad);

end

