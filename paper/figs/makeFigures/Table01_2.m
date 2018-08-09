%% Function to produce Table 01 
%  Table one summarizes the inversions shown in Figure 05
function Table01
ofile = '../Table01';
ifsave = 0;

caption = 'Details of the synthetic tests in Figure 5. Final model parameters for OBSrange inversions are the average of 1000 bootstrap iterations and are omitted if held fixed during the inversion. Parameters $x$ and $y$ are displayed as distance from the drop location.';

%% load 
data_dirs = {
    '2_OUT_wcorr_xrec';
    '1_OUT_nocorr';
    '7_OUT_nocorr_noellipsoid';
    '3_OUT_nocorr_TAT';
    '4_OUT_nocorr_Vp';
    '5_OUT_nocorr_Z';
    '6_OUT_nocorr_TAT_Vp_Z';
    '8_SIO_compare_nobads';
    '9_SIO_compare_wbads';
    };

synth_dirs = {
    '2_OUT_wcorr_xrec';
    '1_OUT_nocorr';
    '7_OUT_wcorr_xrec_noellipsoid';
    '3_OUT_wcorr_xrec_TAT';
    '4_OUT_wcorr_xrec_Vp';
    '5_OUT_wcorr_xrec_Z';
    '6_OUT_wcorr_xrec_TAT_Vp_Z';
    '8_SIO_compare_nobads';
    '9_SIO_compare_wbads';
    };

xlabels = {
    'OBSrange';
    'No Doppler';
    'No Ellipsoid';
    'XYZ$\mathbf{V_p}$';
    'XYZ$\mathbf{\tau}$';
    'XY$\mathbf{\tau V_p}$';
    'XY';
    'SIOgs';
    'SIOgs no QC';
    };

symbols = {
    'pk';
    'ok';
    'ok';
    'ok';
    'ok';
    'ok';
    'ok';
    'pk';
    'ok';
    };

sizes = [
    20;
    14;
    14;
    14;
    14;
    14;
    14;
    20;
    14 ];

%           X Y Z TAT Vp
issolve = [ 1 1 1  1  1;
            1 1 1  1  1;
            1 1 1  1  1;
            1 1 1  0  1;
            1 1 1  1  0;
            1 1 0  1  1;
            1 1 0  0  0;
            1 1 0  0  0;
            1 1 0  0  0 ];
        
%   Doppler Ellipsoid BadPings      
info = {
    'Yes', 'Yes', 'No';
    'No',  'Yes', 'No';
    'Yes', 'No',  'No';
    'Yes', 'Yes', 'No';
    'Yes', 'Yes', 'No';
    'Yes', 'Yes', 'No';
    'Yes', 'Yes', 'No';
    'No' , 'No' , 'No';
    'No' , 'No' , 'Yes';
    };

% swap yes/no of Bad pings (remove bad pings)
for ii = 1:size(info,1)
    if strcmp(info{ii,3},'Yes')
        info{ii,3} = 'No';
    else
        info{ii,3} = 'Yes';
    end
end
        
method = {
    'OBSrange';
    'OBSrange';
    'OBSrange';
    'OBSrange';
    'OBSrange';
    'OBSrange';
    'OBSrange';
    'Grid Search';
    'Grid Search';
    };

synth_path = '../figdata/PacificORCA_synthtest4/OUT_OBSrange';
% Load synthetic
trudata = load('../figdata/PacificORCA_synthtest4/trudata_syn12.mat');

fmt = '%.f';
fmt2 = '%.1f';
for ifil = 1:length(synth_dirs)
    synth_mat = dir(fullfile(synth_path,synth_dirs{ifil},'mats/*.mat'));
    synth = load(fullfile(synth_path,synth_dirs{ifil},'mats',synth_mat.name));
    
    %%%%%% INITIAL
    initial.x_sta{ifil} = num2str(0,fmt);
    initial.y_sta{ifil} = num2str(0,fmt);
    initial.z_sta{ifil} = num2str(synth.datamat.drop_lonlatz(3),fmt);
    if ~(ifil == 8 || ifil == 9)
        initial.TAT{ifil} = num2str(synth.datamat.par.TAT_start*1000,fmt2);
        initial.Vp{ifil} = num2str(synth.datamat.par.vp_w,fmt);
    else
        initial.TAT{ifil} = num2str(0.013,fmt2);
        initial.Vp{ifil} = num2str(1500,fmt);
    end
    
    %%%%%% FINAL
    final.x_sta{ifil} = num2str(median(synth.datamat.x_sta_bs),fmt);
    final.y_sta{ifil} = num2str(median(synth.datamat.y_sta_bs),fmt);
    if ~issolve(ifil,3)
        final.z_sta{ifil} = '-';
    else
        final.z_sta{ifil} = num2str(median(synth.datamat.z_sta_bs),fmt);
    end
    if ~issolve(ifil,4)
        final.TAT{ifil} = '-';
    else
        final.TAT{ifil} = num2str(median(synth.datamat.TAT_bs*1000),fmt2);
    end
    if ~issolve(ifil,5)
        final.Vp{ifil} = '-';
    else
        final.Vp{ifil} = num2str(median(synth.datamat.V_w_bs),fmt);
    end
    
    %%%%%% TRUE
    true.x_sta{ifil} = num2str(trudata.obs_location_xyz(1)*1000,fmt);
    true.y_sta{ifil} = num2str(trudata.obs_location_xyz(2)*1000,fmt);
    true.z_sta{ifil} = num2str(-trudata.obs_location_xyz(3)*1000,fmt);
    true.TAT{ifil} = num2str(trudata.tat*1000,fmt2);
    true.Vp{ifil} = num2str(trudata.vp_actual*1000,fmt);
    
    %%%%%% RMS results
    misfit_xsta_bs(:,ifil) = synth.datamat.x_sta_bs - trudata.obs_location_xyz(1)*1000;
    misfit_ysta_bs(:,ifil) = synth.datamat.y_sta_bs - trudata.obs_location_xyz(2)*1000;
    misfit_zsta_bs(:,ifil) = synth.datamat.z_sta_bs - (-trudata.obs_location_xyz(3)*1000);
    misfit.x_sta{ifil} = num2str(rms(misfit_xsta_bs(:,ifil)),fmt2);
    misfit.y_sta{ifil} = num2str(rms(misfit_ysta_bs(:,ifil)),fmt2);
    misfit.z_sta{ifil} = num2str(rms(misfit_zsta_bs(:,ifil)),fmt2);
    RMS_data(ifil) = mean(synth.datamat.E_rms);
    if ~(ifil == 8 || ifil == 9)
        misfit_TAT_bs(:,ifil) = synth.datamat.TAT_bs - trudata.tat;
        misfit.TAT{ifil} = num2str(rms(misfit_TAT_bs(:,ifil))*1000,fmt2);
        misfit_Vp_bs(:,ifil) = synth.datamat.V_w_bs - trudata.vp_actual*1000;
        misfit.Vp{ifil} = num2str(rms(misfit_Vp_bs(:,ifil)),fmt2);
    else
        misfit.TAT{ifil} = num2str(rms(0.013 - trudata.tat)*1000,fmt2);
        misfit.Vp{ifil} = num2str(rms(1500 - trudata.vp_actual*1000),fmt2);
    end
    misfit_r_xy_bs(:,ifil) = sqrt( misfit_xsta_bs(:,ifil).^2 + misfit_ysta_bs(:,ifil).^2 );
    misfit.r_xy{ifil} = num2str(rms(misfit_r_xy_bs(:,ifil)),fmt2);
   
end


%% ---------------------   Make Latex Table   ---------------------   
fid = fopen([ofile,'.tex'],'w');

% Write Header
fprintf(fid,'\\renewcommand{\\arraystretch}{1.4}'); % row padding
fprintf(fid,'\\begin{table}\n');
fprintf(fid,'\\caption{%s}\n',caption);
fprintf(fid,'\\centering\n');
fprintf(fid,'\\resizebox{\\textwidth}{!}{\n'); % shrink to page width
fprintf(fid,'\\begin{tabular}{c || l c | l c c c c c}\n');
fprintf(fid,'\\textbf{Model Name} &  &  &  & $\\mathbf{x}$ \\textbf{(m)} & $\\mathbf{y}$ \\textbf{(m)} & $\\mathbf{z}$ \\textbf{(m)} & $\\mathbf{\\tau}$ \\textbf{(ms)} & $\\mathbf{V_p}$ \\textbf{(m/s)} \\\\ \n');
fprintf(fid,'\\hline\n');

% Write data
for ifil = 1:length(synth_dirs)
    fprintf(fid,'\\hline\n');
    fprintf(fid,'\\multirow{4}{*}{\\textbf{%s}} & \\textbf{method} & %s & \\textbf{initial} & %s & %s & %s & %s & %s \\\\ \n',xlabels{ifil},method{ifil},initial.x_sta{ifil},initial.y_sta{ifil},initial.z_sta{ifil},initial.TAT{ifil},initial.Vp{ifil});
    fprintf(fid,'\\multirow{4}{*}{} & $\\mathbf{\\delta T}$ \\textbf{correction} & %s & \\textbf{final}& %s & %s & %s & %s & %s \\\\ \n',info{ifil,1},final.x_sta{ifil},final.y_sta{ifil},final.z_sta{ifil},final.TAT{ifil},final.Vp{ifil});
    fprintf(fid,'\\multirow{4}{*}{} & \\textbf{ellipsoid correction} & %s & \\textbf{true}& %s & %s & %s & %s & %s \\\\ \n',info{ifil,2},true.x_sta{ifil},true.y_sta{ifil},true.z_sta{ifil},true.TAT{ifil},true.Vp{ifil});
    fprintf(fid,'\\multirow{4}{*}{} & \\textbf{remove bad data} & %s & \\textbf{RMS} & %s & %s & %s & %s & %s \\\\ \n',info{ifil,3},misfit.x_sta{ifil},misfit.y_sta{ifil},misfit.z_sta{ifil},misfit.TAT{ifil},misfit.Vp{ifil});
end

% Write footer
fprintf(fid,'\\hline\n');
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'}\n'); % shrink to page width
fprintf(fid,'\\label{table:compare_tool}\n');
fprintf(fid,'\\end{table}\n');

fclose(fid);



end

