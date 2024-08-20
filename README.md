# surfviz

MATLAB scripts for visualization of surface-mapped scalar values. Largely inspired by visualization functions in `CAT12` (http://www.neuro.uni-jena.de/cat/) but more useful in serial/parallel batch processing (i.e., no user-interaction required).

## Features
- Various FreeSurfer (http://freesurfer.net/) surface models (white, pial, inflated, semi-inflated, sphere) can be used. 
- OpenGL is required for transparent overlay. Otherwise binary curvature won't be rendered.
- Multiple figures can be created in the background (MATLAB won't bother you during batch processing).
- (NEW) Quick (0.07 sec per `fsaverage` hemisphere) isocontour plotting using an isocontour algorithm from `ft_plot_topo3d()` from FieldTrip (http://www.fieldtriptoolbox.org/) and an edge-connecting algorithm.

## Install
Download/fetch the package and add all subdirectories to MATLAB path:
```Matlab
addpath(genpath(DIRECTORY_WHERE_YOU_SAVED_FILES))
```
Replace "DIRECTORY_WHERE_YOU_SAVED_FILES" with an actual path (e.g., '/Users/me/Documents/MATLAB/surfviz/').

## Demo
```Matlab
subdir = fullfile(getenv('FREESURFER_HOME'),'subjects'); % Find FREESURFER default subject directory
surfs = fsss_read_all_FS_surfs('fsaverage',subdir); % Reading template surfaces & lots of things in the directory
fsss_view(surfs, surfs.THNS) % for example, thickness, areas, white/pial curvature, sulcal depth, annotations, ...
```
![](https://github.com/solleo/surfviz/blob/master/images/demo1.png)

You can also set the surface to visualize, layout, and a threshold via `cfg` structure:
```Matlab
cfg = struct('basesurf','pial', 'layout','1x2', 'thres',3); % just an arbitrary threshold for demo
fsss_view(surfs, surfs.thns, cfg)
```
![](https://github.com/solleo/surfviz/blob/master/images/demo2.2.png)

... and some captions and colorschemes: 
```Matlab
thns_z = {zscore(surfs.thns{1}), zscore(surfs.thns{2})};
cfg = struct('basesurf','inflated', 'layout','2x2', 'thres',1, ...
  'colorbartitle','Rel. Ctx. Thns.', 'colorbarxlabel','Z-score', ...
  'colorscheme','yellowblue');
fsss_view(surfs, thns_z, cfg)
```
![](https://github.com/solleo/surfviz/blob/master/images/demo4.3.png)

... finally you can create as many figures as you want without Matalb taking away your attention (figure's visibility will be `off`):
```Matlab
fsss_view(surfs, surfs.THNS, struct('demo.png')) % Nothing pops up but it creates a PNG file _silently_
ls('demo.png') % check the file was created!
```

(NEW) With an edge-connecting algorithm, isocontours can be plotted very quickly (0.07 sec per 160k-vert surface). Isocurvature can be very useful for unthresholded maps to mark the boundary between gyri and sulci. Any arbitrary contours (significant clusters, manual/atals-based ROIs) can be also overlaid:
```Matlab
subdir = fullfile(getenv('FREESURFER_HOME'),'subjects'); % Find FREESURFER default subject directory
surfs = fsss_read_all_FS_surfs('fsaverage',subdir, struct('isocurv','inflated')); % Computes isocurvature line groups when loading
fsss_view(surfs, surfs.thns)
```
![](https://github.com/solleo/surfviz/blob/master/images/demo5.png)

Or perhaps the HCP-MMP atlas ([projected on fsaverage](https://figshare.com/articles/HCP-MMP1_0_projected_on_fsaverage/3498446))?
```Matlab
fsavg = fsss_read_all_FS_surfs('fsaverage');
labels = {}; cots = {};
[~, labels{1}, cots{1}] = read_annotation('lh.HCP-MMP1.annot', 0);
[~, labels{2}, cots{2}] = read_annotation('rh.HCP-MMP1.annot', 0);
[fsavg,cots] = fsss_isoclus(fsavg, labels, struct('cots',{cots})); % this takes a while for granular annotations
fsss_view(fsavg, fsavg.thns, struct('colormap',gray, 'isocluslinewidth',1,'isocluscolors',cots)
```
![](https://github.com/solleo/surfviz/blob/master/images/demo6.png)

See documentation for more information:
```Matlab
doc fsss_view
```

## Publications
The current or previous version of this repo has been used in these publications:
- [Kim et al., 2024, Cerebral Cortex](https://doi.org/10.1093/cercor/bhae155)
- [Weidacker, Kim et al., 2020, Psychological Medicine](https://doi.rog/10.1017/S0033291720003852)
- [Kim, Overath et al., 2022, NeuroImage](https://doi.org/10.1016/j.neuroimage.2022.118879)
- [Nord, Kim et al., 2019, Neuropsychopharmacology](https://doi.org/10.1038/s41386-019-0343-6)
- [Zhang et al., 2019, Brain Stimulation](https://doi.org/10.1016/j.brs.2019.05.010)

(cc) Seung-Goo Kim, 2019-2024.
