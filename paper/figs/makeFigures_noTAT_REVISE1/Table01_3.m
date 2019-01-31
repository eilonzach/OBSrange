%% Function to produce Table 01 
%  Table one summarizes the inversions shown in Figure 05
function Table01
ofile = '../Table01';
ifsave = 0;

caption = 'Details of the synthetic tests in Figure 3 for a \textit{PACMAN} survey of radius 1~Nm and 5050~m instrument depth. Final model parameters for OBSrange inversions are the average of 1000 bootstrap iterations. Parameters that are held fixed during the inversion are denoted in italics and their final values omitted. Parameters $x$ and $y$ are displayed as distance from the drop location.';

%% load 
iSIO = [7 8]+1;
data_dirs = {
    '2_OUT_wcorr_xrec';
    '1_OUT_nocorr';
    '7_OUT_nocorr_noellipsoid';
%     '3_OUT_nocorr_TAT';
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
    '10_OUT_wcorr_xrec_noGPScorr';
%     '3_OUT_wcorr_xrec_TAT';
    '4_OUT_wcorr_xrec_Vp';
    '5_OUT_wcorr_xrec_Z';
    '6_OUT_wcorr_xrec_TAT_Vp_Z';
    '8_SIO_compare_nobads';
    '9_SIO_compare_wbads';
    };

% xlabels = {
%     'OBSrange';
%     'No Doppler';
%     'No Ellipsoid';
%     'XYZ$V_p$';
%     'XYZ$\tau$';
%     'XY$\tau V_p$';
%     'XY';
%     'SIOgs';
%     'SIOgs no QC';
%     };

xlabels = {
    '(1) OBSrange';
    '(2) No Doppler';
    '(3) No Ellipsoid';
    '(4) No GPS';
%     'Fix-$\tau$';
    '(5) Fix-$V_p$';
    '(6) Fix-Z';
    '(7) XY-only';
    '(8) SIOgs';
    '(9) SIOgs no QC';
    };

symbols = {
    'pk';
    'ok';
    'ok';
    'ok';
%     'ok';
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
%     14;
    14;
    14;
    14;
    20;
    14 ];

%           X Y Z TAT Vp
issolve = [ 1 1 1  1  1;
            1 1 1  1  1;
            1 1 1  1  1;
            1 1 1  1  1;
%             1 1 1  0  1;
            1 1 1  1  0;
            1 1 0  1  1;
            1 1 0  0  0;
            1 1 0  0  0;
            1 1 0  0  0 ];
        
%   Doppler Ellipsoid BadPings GPS
info = {
    'Yes', 'Yes', 'No', 'Yes';
    'No',  'Yes', 'No', 'Yes';
    'Yes', 'No',  'No', 'Yes';
    'Yes', 'Yes', 'No', 'No';
%     'Yes', 'Yes', 'No';
    'Yes', 'Yes', 'No', 'Yes';
    'Yes', 'Yes', 'No', 'Yes';
    'Yes', 'Yes', 'No', 'Yes';
    'No' , 'No' , 'No', 'No';
    'No' , 'No' , 'Yes','No';
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
%     'OBSrange';
    'OBSrange';
    'OBSrange';
    'OBSrange';
    'Grid Search';
    'Grid Search';
    };

% synth_path = '../figdata/PacificORCA_synthtest4_noTAT/OUT_OBSrange';
% % Load synthetic
% trudata = load('../figdata/PacificORCA_synthtest4_noTAT/trudata_syn12.mat');

% REVISION1
synth_path = '../figdata/PacificORCA_synthtest4_REVISION1_GPScorr/OUT_OBSrange';
trudata = load('../figdata/PacificORCA_synthtest4_REVISION1_GPScorr/trudata_syn12_z5000m_fr10.mat');

fmt = '%.f';
fmt2 = '%.1f';
for ifil = 1:length(synth_dirs)
    synth_mat = dir(fullfile(synth_path,synth_dirs{ifil},'mats/*.mat'));
    synth = load(fullfile(synth_path,synth_dirs{ifil},'mats',synth_mat.name));
    
    %%%%%% INITIAL
    initial.x_sta{ifil} = num2str(0,fmt);
    initial.y_sta{ifil} = num2str(0,fmt);
    initial.z_sta{ifil} = num2str(synth.datamat.drop_lonlatz(3),fmt);
    if ~(ifil == iSIO(1) || ifil == iSIO(2))
%         initial.TAT{ifil} = num2str(synth.datamat.par.TAT_start*1000,fmt2);
        initial.Vp{ifil} = num2str(synth.datamat.par.vp_w,fmt);
    else
%         initial.TAT{ifil} = num2str(0.013*1000,fmt2);
        initial.Vp{ifil} = num2str(1500,fmt);
    end
    
    %%%%%% FINAL
    final.x_sta{ifil} = num2str(median(synth.datamat.x_sta_bs),fmt);
    final.y_sta{ifil} = num2str(median(synth.datamat.y_sta_bs),fmt);
    if ~issolve(ifil,3)
        final.z_sta{ifil} = '-';
        initial.z_sta{ifil} = ['$\mathit{',initial.z_sta{ifil},'}$'];
    else
        final.z_sta{ifil} = num2str(median(synth.datamat.z_sta_bs),fmt);
    end
    if ~issolve(ifil,4)
%         final.TAT{ifil} = '-';
%         initial.TAT{ifil} = ['$\mathit{',initial.TAT{ifil},'}$'];
    else
%         final.TAT{ifil} = num2str(median(synth.datamat.TAT_bs*1000),fmt2);
    end
    if ~issolve(ifil,5)
        final.Vp{ifil} = '-';
        initial.Vp{ifil} = ['$\mathit{',initial.Vp{ifil},'}$'];
    else
        final.Vp{ifil} = num2str(median(synth.datamat.V_w_bs),fmt);
    end
    
    %%%%%% TRUE
    true.x_sta{ifil} = num2str(trudata.obs_location_xyz(1)*1000,fmt);
    true.y_sta{ifil} = num2str(trudata.obs_location_xyz(2)*1000,fmt);
    true.z_sta{ifil} = num2str(-trudata.obs_location_xyz(3)*1000,fmt);
%     true.TAT{ifil} = num2str(trudata.tat*1000,fmt2);
    true.Vp{ifil} = num2str(trudata.vp_actual*1000,fmt);
    
    %%%%%% RMS results
    misfit_xsta_bs(:,ifil) = synth.datamat.x_sta_bs - trudata.obs_location_xyz(1)*1000;
    misfit_ysta_bs(:,ifil) = synth.datamat.y_sta_bs - trudata.obs_location_xyz(2)*1000;
    misfit_zsta_bs(:,ifil) = synth.datamat.z_sta_bs - (-trudata.obs_location_xyz(3)*1000);
    misfit.x_sta{ifil} = num2str(rms(misfit_xsta_bs(:,ifil)),fmt2);
    misfit.y_sta{ifil} = num2str(rms(misfit_ysta_bs(:,ifil)),fmt2);
    misfit.z_sta{ifil} = num2str(rms(misfit_zsta_bs(:,ifil)),fmt2);
    RMS_data(ifil) = mean(synth.datamat.E_rms);
    if ~(ifil == iSIO(1) || ifil == iSIO(2))
%         misfit_TAT_bs(:,ifil) = synth.datamat.TAT_bs - trudata.tat;
%         misfit.TAT{ifil} = num2str(rms(misfit_TAT_bs(:,ifil))*1000,fmt2);
        misfit_Vp_bs(:,ifil) = synth.datamat.V_w_bs - trudata.vp_actual*1000;
        misfit.Vp{ifil} = num2str(rms(misfit_Vp_bs(:,ifil)),fmt2);
    else
%         misfit.TAT{ifil} = num2str(rms(0.013 - trudata.tat)*1000,fmt2);
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
fprintf(fid,'\\begin{tabular}{l | c c c c c | l c c c c}\n');
fprintf(fid,'                     &                  & \\textbf{Doppler}     & \\textbf{ellipsoid}  & \\textbf{GPS}         & \\textbf{remove}     & & $\\mathbf{x}$ & $\\mathbf{y}$ & $\\mathbf{z}$ & $\\mathbf{V_p}$ \\\\ \n');
fprintf(fid,'\\textbf{Model Name} & \\textbf{method} & \\textbf{correction}  & \\textbf{correction} & \\textbf{correction}  & \\textbf{bad data}   & & \\textbf{(m)} & \\textbf{(m)} & \\textbf{(m)} & \\textbf{(m/s)} \\\\ \n');
fprintf(fid,'\\hline\n');

% Write data
for ifil = 1:length(synth_dirs)
    fprintf(fid,'\\hline\n');
    fprintf(fid,'\\multirow{4}{*}{\\textbf{%s}} & \\multirow{4}{*}{%s} & \\multirow{4}{*}{%s} & \\multirow{4}{*}{%s} & \\multirow{4}{*}{%s} & \\multirow{4}{*}{%s} & \\textbf{initial} & %s & %s & %s & %s \\\\ \n',xlabels{ifil},method{ifil},info{ifil,1}, info{ifil,2}, info{ifil,4},info{ifil,3},initial.x_sta{ifil},initial.y_sta{ifil},initial.z_sta{ifil},initial.Vp{ifil});
    fprintf(fid,'\\multirow{4}{*}{} & \\multirow{4}{*}{} & \\multirow{4}{*}{} & \\multirow{4}{*}{} & \\multirow{4}{*}{} & \\multirow{4}{*}{} & \\textbf{final}& %s & %s & %s & %s \\\\ \n',final.x_sta{ifil},final.y_sta{ifil},final.z_sta{ifil},final.Vp{ifil});
    fprintf(fid,'\\multirow{4}{*}{} & \\multirow{4}{*}{} & \\multirow{4}{*}{} & \\multirow{4}{*}{} & \\multirow{4}{*}{} & \\multirow{4}{*}{} & \\textbf{true}& %s & %s & %s & %s \\\\ \n',true.x_sta{ifil},true.y_sta{ifil},true.z_sta{ifil},true.Vp{ifil});
    fprintf(fid,'\\multirow{4}{*}{} & \\multirow{4}{*}{} & \\multirow{4}{*}{} & \\multirow{4}{*}{} & \\multirow{4}{*}{} & \\multirow{4}{*}{} & \\textbf{RMS} & %s & %s & %s & %s \\\\ \n',misfit.x_sta{ifil},misfit.y_sta{ifil},misfit.z_sta{ifil},misfit.Vp{ifil});
end

% Write footer
fprintf(fid,'\\hline\n');
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'}\n'); % shrink to page width
fprintf(fid,'\\label{table:compare_tool}\n');
fprintf(fid,'\\end{table}\n');

fclose(fid);



end

