%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start-to-end MATLAB implementation of the protocol
%
% Ensure FieldTrip is correcty added to the MATLAB path:
% addpath <path to fieldtrip home directory>
% ft_defaults
%
% The data is available at: ftp://ftp.fieldtriptoolbox.org/pub/ ...
% fieldtrip/tutorial/SubjectUCI29.zip
%
% This script is part of Stolk, Griffin et al., Integrated analysis
% of anatomical and electrophysiological human intracranial data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath C:\Users\cours\Downloads\fieldtrip-master
ft_defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subjID = 'SubjectPD22N008';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preprocessing of the anatomical MRI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mri = ft_read_mri([subjID '_MR_acpc.nii']); % we used the dcm series
ft_determine_coordsys(mri);
cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'acpc';
mri_acpc = ft_volumerealign(cfg, mri);
cfg = [];
cfg.filename = [subjID '_MR_acpc'];
cfg.filetype = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, mri_acpc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FreeSurfer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fshome = "C:\Users\cours\Downloads\freesurfer\";
subdir = pwd; % present working directory
command = sprintf(['export FREESURFER_HOME=%s; ' ...
    'source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' ...
    'mri_convert -c -oc 0 0 0 %s %s; ' ...
    'recon-all -i %s -s freesurfer -sd %s -all'], ...
    fshome, mrfile, fullfile(subdir, 'tmp.nii'), fullfile(subdir, 'tmp.nii'), subdir);
system(command);

% mrfile = [subdir filesep subjID '_MR_acpc.nii'];
% system(['export FREESURFER_HOME=' fshome '; ' ...
%  'source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' ...
%  'mri_convert -c -oc 0 0 0 ' mrfile ' ' [subdir '/tmp.nii'] '; ' ...
%  'recon-all -i ' [subdir '/tmp.nii'] ' -s ' 'freesurfer' ' -sd ' ...
%  subdir ' -all'])
fsmri_acpc = ft_read_mri('freesurfer/mri/T1.mgz');
fsmri_acpc.coordsys = 'acpc';
pial_lh = ft_read_headshape('freesurfer/surf/lh.pial.T1');
pial_lh.coordsys = 'acpc';
%1;
ft_plot_mesh(pial_lh);
material dull; lighting gouraud; camlight;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preprocessing of the anatomical CT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ct = ft_read_mri([subjID '_CT_acpc_f.nii']); % we used the dcm series
cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'ctf';
ct_ctf = ft_volumerealign(cfg, ct);
ct_acpc = ft_convert_coordsys(ct_ctf, 'acpc');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fusion of the CT with the MRI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg = [];
cfg.method = 'spm';
cfg.spmversion = 'spm12';
cfg.coordsys = 'acpc';
cfg.viewresult = 'yes';
ct_acpc_f = ft_volumerealign(cfg, ct_acpc, fsmri_acpc);
cfg = [];
cfg.filename = [subjID '_CT_acpc_f'];
cfg.filetype = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, ct_acpc_f);
print([subjID '_CT_acpc_f.png'], '-dpng');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% electrode placement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([subjID '_hdr.mat']); % we used ft_read_header
cfg = [];
cfg.channel = hdr.label;
elec_acpc_f = ft_electrodeplacement(cfg, ct_acpc_f, fsmri_acpc);
ft_plot_ortho(fsmri_acpc.anatomy, 'transform', ...
 fsmri_acpc.transform, 'style', 'intersect');
ft_plot_sens(elec_acpc_f, 'label', 'on', 'fontcolor', 'w')
print([subjID '_elec_acpc_f.png'], '-dpng');
save([subjID '_elec_acpc_f.mat'], 'elec_acpc_f');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% brain shift compensation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2
cfg = [];
cfg.method = 'cortexhull';
cfg.headshape = 'freesurfer/surf/lh.pial';
cfg.fshome = '/Applications/freesurfer';
hull_lh = ft_prepare_mesh(cfg);
save([subjID '_hull_lh.mat'], 'mesh');
elec_acpc_fr = elec_acpc_f;
grids = {'LPG*', 'LTG*'};
for g = 1:numel(grids)
 cfg = [];
 cfg.channel = grids{g};
 cfg.keepchannel = 'yes';
 cfg.elec = elec_acpc_fr;
 cfg.method = 'headshape';
 cfg.headshape = hull_lh;
 cfg.warp = 'dykstra2012';
 cfg.feedback = 'yes';
 elec_acpc_fr = ft_electroderealign(cfg);
end
ft_plot_mesh(pial_lh);
ft_plot_sens(elec_acpc_fr);
view([-55 10]); material dull; lighting gouraud; camlight
save([subjID '_elec_acpc_fr.mat'], 'elec_acpc_fr');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% volume-based registration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ftver, ftpath] = ft_version;
cfg = [];
cfg.nonlinear = 'yes';
cfg.spmversion = 'spm12';
fsmri_mni = ft_volumenormalise(cfg, fsmri_acpc);
elec_mni_frv = elec_acpc_fr;
elec_mni_frv.elecpos = ft_warp_apply(fsmri_mni.params, ...
 elec_acpc_fr.elecpos, 'individual2sn');
elec_mni_frv.chanpos = ft_warp_apply(fsmri_mni.params, ...
 elec_acpc_fr.chanpos, 'individual2sn');
elec_mni_frv.coordsys = 'mni';
save([subjID '_elec_mni_frv.mat'], 'elec_mni_frv');
load([ftpath filesep 'template/anatomy/surface_pial_left.mat']);
ft_plot_mesh(mesh);
ft_plot_sens(elec_mni_frv);
view([-90 20]); material dull; lighting gouraud; camlight;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% surface-based registration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
3
cfg = [];
cfg.channel = {'LPG*', 'LTG*'};
cfg.elec = elec_acpc_fr;
cfg.method = 'headshape';
cfg.headshape = 'freesurfer/surf/lh.pial';
cfg.warp = 'fsaverage';
cfg.fshome = '/Applications/freesurfer';
elec_fsavg_frs = ft_electroderealign(cfg);
save([subjID '_elec_fsavg_frs.mat'], 'elec_fsavg_frs');
fspial_lh = ft_read_headshape( ...
 '/Applications/freesurfer/subjects/fsaverage/surf/lh.pial');
fspial_lh.coordsys = 'fsaverage';
ft_plot_mesh(fspial_lh);
ft_plot_sens(elec_fsavg_frs);
view([-90 20]); material dull; lighting gouraud; camlight;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% anatomical labeling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
atlas = ft_read_atlas([ftpath filesep ...
 'template/atlas/aal/ROI_MNI_V4.nii']);
cfg = [];
cfg.roi = elec_mni_frv.chanpos( ...
 match_str(elec_mni_frv.label, 'LHH1'),:);
cfg.atlas = atlas;
cfg.inputcoord = 'mni';
cfg.output = 'label';
labels = ft_volumelookup(cfg, atlas);
[~, indx] = max(labels.count);
labels.name(indx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data inspection and artifact rejection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([subjID '_data.mat'], 'data'); % we used ft_preprocessing
cfg = [];
cfg.viewmode = 'vertical';
cfg = ft_databrowser(cfg, data);
data = ft_rejectartifact(cfg, data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% re-referencing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg = [];
cfg.channel = {'LPG*', 'LTG*'};
cfg.reref = 'yes';
cfg.refchannel = 'all';
cfg.refmethod = 'avg';
4
reref_grids = ft_preprocessing(cfg, data);
depths = {'RAM*', 'RHH*', 'RTH*', 'ROC*', 'LAM*', 'LHH*', 'LTH*'};
for d = 1:numel(depths)
 cfg = [];
 cfg.channel = ft_channelselection(depths{d}, data.label);
 cfg.reref = 'yes';
 cfg.refchannel = 'all';
 cfg.refmethod = 'bipolar';
 cfg.updatesens = 'yes';
 reref_depths{d} = ft_preprocessing(cfg, data);
end
cfg = [];
cfg.appendsens = 'yes';
reref = ft_appenddata(cfg, reref_grids, reref_depths{:});
save([subjID '_reref.mat'], 'reref');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time-frequency analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg = [];
cfg.method = 'mtmconvol';
cfg.foi = 5:5:200;
cfg.toi = -.3:0.01:.8;
cfg.t_ftimwin = ones(length(cfg.foi),1).*0.2;
cfg.taper = 'hanning';
cfg.output = 'pow';
cfg.keeptrials = 'no';
freq = ft_freqanalysis(cfg, reref);
save([subjID '_freq.mat'], 'freq');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interactive plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg = [];
cfg.headshape = pial_lh;
cfg.projection = 'orthographic';
cfg.channel = {'LPG*', 'LTG*'};
cfg.viewpoint = 'left';
cfg.mask = 'convex';
cfg.boxchannel = {'LTG30', 'LTG31'};
lay = ft_prepare_layout(cfg, freq);
cfg = [];
cfg.baseline = [-.3 -.1];
cfg.baselinetype = 'relchange';
freq_blc = ft_freqbaseline(cfg, freq);
cfg = [];
cfg.layout = lay;
cfg.showoutline = 'yes';
5
ft_multiplotTFR(cfg, freq_blc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECoG data representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg = [];
cfg.frequency = [70 150];
cfg.avgoverfreq = 'yes';
cfg.latency = [0 0.8];
cfg.avgovertime = 'yes';
freq_sel = ft_selectdata(cfg, freq_blc);
cfg = [];
cfg.funparameter = 'powspctrm';
cfg.funcolorlim = [-.5 .5];
cfg.method = 'surface';
cfg.interpmethod = 'sphere_weighteddistance';
cfg.sphereradius = 8;
cfg.camlight = 'no';
ft_sourceplot(cfg, freq_sel, pial_lh);
view([-90 20]); material dull; lighting gouraud; camlight;
ft_plot_sens(elec_acpc_fr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEEG data representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
atlas = ft_read_atlas('freesurfer/mri/aparc+aseg.mgz');
atlas.coordsys = 'acpc';
cfg = [];
cfg.inputcoord = 'acpc';
cfg.atlas = atlas;
cfg.roi = {'Right-Hippocampus', 'Right-Amygdala'};
mask_rha = ft_volumelookup(cfg, atlas);
seg = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
seg.brain = mask_rha;
cfg = [];
cfg.method = 'iso2mesh';
cfg.numvertices = 10000;
cfg.radbound = 2;
cfg.maxsurf = 0;
cfg.tissue = 'brain';
cfg.smooth = 3;
cfg.spmversion = 'spm12';
mesh_rha = ft_prepare_mesh(cfg, seg);
cfg = [];
cfg.channel = {'RAM*', 'RTH*', 'RHH*'};
freq_sel2 = ft_selectdata(cfg, freq_sel);
cfg = [];
cfg.funparameter = 'powspctrm';
cfg.funcolorlim = [-.5 .5];
6
cfg.method = 'cloud';
cfg.slice = '3d';
cfg.nslices = 2;
cfg.facealpha = .25;
ft_sourceplot(cfg, freq_sel2, mesh_rha);
view([120 40]); lighting gouraud; camlight;
cfg.slice = '2d';
ft_sourceplot(cfg, freq_sel2, mesh_rha);
Published with MATLAB® R2017b
7%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start-to-end MATLAB implementation of the protocol
%
% Ensure FieldTrip is correcty added to the MATLAB path:
% addpath <path to fieldtrip home directory>
% ft_defaults
%
% The data is available at: ftp://ftp.fieldtriptoolbox.org/pub/ ...
% fieldtrip/tutorial/SubjectUCI29.zip
%
% This script is part of Stolk, Griffin et al., Integrated analysis
% of anatomical and electrophysiological human intracranial data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subjID = 'SubjectUCI29';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preprocessing of the anatomical MRI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mri = ft_read_mri([subjID '_MR_acpc.nii']); % we used the dcm series
ft_determine_coordsys(mri);
cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'acpc';
mri_acpc = ft_volumerealign(cfg, mri);
cfg = [];
cfg.filename = [subjID '_MR_acpc'];
cfg.filetype = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, mri_acpc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FreeSurfer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fshome = '/Applications/freesurfer';
subdir = pwd; % present working directory
mrfile = [subdir filesep subjID '_MR_acpc.nii'];
system(['export FREESURFER_HOME=' fshome '; ' ...
 'source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' ...
 'mri_convert -c -oc 0 0 0 ' mrfile ' ' [subdir '/tmp.nii'] '; ' ...
 'recon-all -i ' [subdir '/tmp.nii'] ' -s ' 'freesurfer' ' -sd ' ...
 subdir ' -all'])
fsmri_acpc = ft_read_mri('freesurfer/mri/T1.mgz');
fsmri_acpc.coordsys = 'acpc';
pial_lh = ft_read_headshape('freesurfer/surf/lh.pial');
pial_lh.coordsys = 'acpc';
1
ft_plot_mesh(pial_lh);
material dull; lighting gouraud; camlight;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preprocessing of the anatomical CT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ct = ft_read_mri([subjID '_CT_acpc_f.nii']); % we used the dcm series
cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'ctf';
ct_ctf = ft_volumerealign(cfg, ct);
ct_acpc = ft_convert_coordsys(ct_ctf, 'acpc');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fusion of the CT with the MRI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg = [];
cfg.method = 'spm';
cfg.spmversion = 'spm12';
cfg.coordsys = 'acpc';
cfg.viewresult = 'yes';
ct_acpc_f = ft_volumerealign(cfg, ct_acpc, fsmri_acpc);
cfg = [];
cfg.filename = [subjID '_CT_acpc_f'];
cfg.filetype = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, ct_acpc_f);
print([subjID '_CT_acpc_f.png'], '-dpng');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% electrode placement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([subjID '_hdr.mat']); % we used ft_read_header
cfg = [];
cfg.channel = hdr.label;
elec_acpc_f = ft_electrodeplacement(cfg, ct_acpc_f, fsmri_acpc);
ft_plot_ortho(fsmri_acpc.anatomy, 'transform', ...
 fsmri_acpc.transform, 'style', 'intersect');
ft_plot_sens(elec_acpc_f, 'label', 'on', 'fontcolor', 'w')
print([subjID '_elec_acpc_f.png'], '-dpng');
save([subjID '_elec_acpc_f.mat'], 'elec_acpc_f');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% brain shift compensation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
2
cfg = [];
cfg.method = 'cortexhull';
cfg.headshape = 'freesurfer/surf/lh.pial';
cfg.fshome = '/Applications/freesurfer';
hull_lh = ft_prepare_mesh(cfg);
save([subjID '_hull_lh.mat'], 'mesh');
elec_acpc_fr = elec_acpc_f;
grids = {'LPG*', 'LTG*'};
for g = 1:numel(grids)
 cfg = [];
 cfg.channel = grids{g};
 cfg.keepchannel = 'yes';
 cfg.elec = elec_acpc_fr;
 cfg.method = 'headshape';
 cfg.headshape = hull_lh;
 cfg.warp = 'dykstra2012';
 cfg.feedback = 'yes';
 elec_acpc_fr = ft_electroderealign(cfg);
end
ft_plot_mesh(pial_lh);
ft_plot_sens(elec_acpc_fr);
view([-55 10]); material dull; lighting gouraud; camlight
save([subjID '_elec_acpc_fr.mat'], 'elec_acpc_fr');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% volume-based registration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ftver, ftpath] = ft_version;
cfg = [];
cfg.nonlinear = 'yes';
cfg.spmversion = 'spm12';
fsmri_mni = ft_volumenormalise(cfg, fsmri_acpc);
elec_mni_frv = elec_acpc_fr;
elec_mni_frv.elecpos = ft_warp_apply(fsmri_mni.params, ...
 elec_acpc_fr.elecpos, 'individual2sn');
elec_mni_frv.chanpos = ft_warp_apply(fsmri_mni.params, ...
 elec_acpc_fr.chanpos, 'individual2sn');
elec_mni_frv.coordsys = 'mni';
save([subjID '_elec_mni_frv.mat'], 'elec_mni_frv');
load([ftpath filesep 'template/anatomy/surface_pial_left.mat']);
ft_plot_mesh(mesh);
ft_plot_sens(elec_mni_frv);
view([-90 20]); material dull; lighting gouraud; camlight;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% surface-based registration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
3
cfg = [];
cfg.channel = {'LPG*', 'LTG*'};
cfg.elec = elec_acpc_fr;
cfg.method = 'headshape';
cfg.headshape = 'freesurfer/surf/lh.pial';
cfg.warp = 'fsaverage';
cfg.fshome = '/Applications/freesurfer';
elec_fsavg_frs = ft_electroderealign(cfg);
save([subjID '_elec_fsavg_frs.mat'], 'elec_fsavg_frs');
fspial_lh = ft_read_headshape( ...
 '/Applications/freesurfer/subjects/fsaverage/surf/lh.pial');
fspial_lh.coordsys = 'fsaverage';
ft_plot_mesh(fspial_lh);
ft_plot_sens(elec_fsavg_frs);
view([-90 20]); material dull; lighting gouraud; camlight;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% anatomical labeling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
atlas = ft_read_atlas([ftpath filesep ...
 'template/atlas/aal/ROI_MNI_V4.nii']);
cfg = [];
cfg.roi = elec_mni_frv.chanpos( ...
 match_str(elec_mni_frv.label, 'LHH1'),:);
cfg.atlas = atlas;
cfg.inputcoord = 'mni';
cfg.output = 'label';
labels = ft_volumelookup(cfg, atlas);
[~, indx] = max(labels.count);
labels.name(indx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data inspection and artifact rejection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([subjID '_data.mat'], 'data'); % we used ft_preprocessing
cfg = [];
cfg.viewmode = 'vertical';
cfg = ft_databrowser(cfg, data);
data = ft_rejectartifact(cfg, data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% re-referencing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg = [];
cfg.channel = {'LPG*', 'LTG*'};
cfg.reref = 'yes';
cfg.refchannel = 'all';
cfg.refmethod = 'avg';
4
reref_grids = ft_preprocessing(cfg, data);
depths = {'RAM*', 'RHH*', 'RTH*', 'ROC*', 'LAM*', 'LHH*', 'LTH*'};
for d = 1:numel(depths)
 cfg = [];
 cfg.channel = ft_channelselection(depths{d}, data.label);
 cfg.reref = 'yes';
 cfg.refchannel = 'all';
 cfg.refmethod = 'bipolar';
 cfg.updatesens = 'yes';
 reref_depths{d} = ft_preprocessing(cfg, data);
end
cfg = [];
cfg.appendsens = 'yes';
reref = ft_appenddata(cfg, reref_grids, reref_depths{:});
save([subjID '_reref.mat'], 'reref');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time-frequency analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg = [];
cfg.method = 'mtmconvol';
cfg.foi = 5:5:200;
cfg.toi = -.3:0.01:.8;
cfg.t_ftimwin = ones(length(cfg.foi),1).*0.2;
cfg.taper = 'hanning';
cfg.output = 'pow';
cfg.keeptrials = 'no';
freq = ft_freqanalysis(cfg, reref);
save([subjID '_freq.mat'], 'freq');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interactive plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg = [];
cfg.headshape = pial_lh;
cfg.projection = 'orthographic';
cfg.channel = {'LPG*', 'LTG*'};
cfg.viewpoint = 'left';
cfg.mask = 'convex';
cfg.boxchannel = {'LTG30', 'LTG31'};
lay = ft_prepare_layout(cfg, freq);
cfg = [];
cfg.baseline = [-.3 -.1];
cfg.baselinetype = 'relchange';
freq_blc = ft_freqbaseline(cfg, freq);
cfg = [];
cfg.layout = lay;
cfg.showoutline = 'yes';
5
ft_multiplotTFR(cfg, freq_blc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECoG data representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg = [];
cfg.frequency = [70 150];
cfg.avgoverfreq = 'yes';
cfg.latency = [0 0.8];
cfg.avgovertime = 'yes';
freq_sel = ft_selectdata(cfg, freq_blc);
cfg = [];
cfg.funparameter = 'powspctrm';
cfg.funcolorlim = [-.5 .5];
cfg.method = 'surface';
cfg.interpmethod = 'sphere_weighteddistance';
cfg.sphereradius = 8;
cfg.camlight = 'no';
ft_sourceplot(cfg, freq_sel, pial_lh);
view([-90 20]); material dull; lighting gouraud; camlight;
ft_plot_sens(elec_acpc_fr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEEG data representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
atlas = ft_read_atlas('freesurfer/mri/aparc+aseg.mgz');
atlas.coordsys = 'acpc';
cfg = [];
cfg.inputcoord = 'acpc';
cfg.atlas = atlas;
cfg.roi = {'Right-Hippocampus', 'Right-Amygdala'};
mask_rha = ft_volumelookup(cfg, atlas);
seg = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
seg.brain = mask_rha;
cfg = [];
cfg.method = 'iso2mesh';
cfg.numvertices = 10000;
cfg.radbound = 2;
cfg.maxsurf = 0;
cfg.tissue = 'brain';
cfg.smooth = 3;
cfg.spmversion = 'spm12';
mesh_rha = ft_prepare_mesh(cfg, seg);
cfg = [];
cfg.channel = {'RAM*', 'RTH*', 'RHH*'};
freq_sel2 = ft_selectdata(cfg, freq_sel);
cfg = [];
cfg.funparameter = 'powspctrm';
cfg.funcolorlim = [-.5 .5];
6
cfg.method = 'cloud';
cfg.slice = '3d';
cfg.nslices = 2;
cfg.facealpha = .25;
ft_sourceplot(cfg, freq_sel2, mesh_rha);
view([120 40]); lighting gouraud; camlight;
cfg.slice = '2d';
ft_sourceplot(cfg, freq_sel2, mesh_rha);
Published with MATLAB® R2017b
7