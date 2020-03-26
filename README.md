# surfviz

MATLAB codes for visualization of surface-mapped scalar values. Largely inspired by visualization functions in `CAT12` (http://www.neuro.uni-jena.de/cat/) but more useful in serial/parallel batch processing (i.e., no user-interaction required).

## Features
- Various FreeSurfer (http://freesurfer.net/) surface models (white, pial, inflated, semi-inflated, sphere) can be used. 
- OpenGL is required for transparent overlay. Otherwise binary curvature won't be rendered.
- Multiple figures can be created in the background (MATLAB won't bother you during batch processing).
- (NEW) Quick (0.07 sec per `fsaverage` hemisphere) isocontour plotting using a modification of `ft_plot_topo3d()` from FieldTrip (http://www.fieldtriptoolbox.org/).

## Install
Download/fetch the package and add all subdirectories to MATLAB path:
```Matlab
addpath(genpath(DIRECTORY_WHERE_YOU_SAVED_FILES))
```

## Demo
```Matlab
subdir = fullfile(getenv('FREESURFER_HOME'),'subjects'); % Find FREESURFER default subject directory
surfs = fsss_read_all_FS_surfs('fsaverage',subdir); % Reading template surfaces & lots of things in the directory
fsss_view(surfs, surfs.THNS) % for example, thickness, areas, white/pial curvature, sulcal depth, annotations, ...
```
![](https://github.com/solleo/surfviz/blob/master/images/demo1.png)

You can also set the surface to visualize, layout, and a threshold via `cfg` structure:
```Matlab
cfg = struct('basesurf','PIAL', 'layout','1x2', 'thres',3); % just an arbitrary threshold for demo
fsss_view(surfs, surfs.THNS, cfg)
```
![](https://github.com/solleo/surfviz/blob/master/images/demo2.2.png)

... and some captions and colorschemes: 
```Matlab
thns_z = {zscore(surfs.THNS{1}), zscore(surfs.THNS{2})};
cfg = struct('basesurf','INFL', 'layout','2x2', 'thres',1, ...
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
surfs = fsss_read_all_FS_surfs('fsaverage',subdir, struct('isocurv','INFL')); % Computes isocurvature line groups when loading
fsss_view(surfs, surfs.THNS)
```
![](https://github.com/solleo/surfviz/blob/master/images/demo5.png)

Or perhaps the HCP-MMP atlas ([projected on fsaverage](https://figshare.com/articles/HCP-MMP1_0_projected_on_fsaverage/3498446))?
```Matlab
fsavg = fsss_read_all_FS_surfs('fsaverage');
labels = {};
[~, labels{1}, ~] = read_annotation( 'lh.HCP-MMP1.annot', 0 );
[~, labels{2}, ~] = read_annotation( 'rh.HCP-MMP1.annot', 0 );
fsavg = fsss_isoclus(fsavg, labels); % this takes a while for granular annotations
```
![](https://github.com/solleo/surfviz/blob/master/images/demo6.png)

See documentation for more information:
```Matlab
doc fsss_view
```

(cc) Seung-Goo Kim, 2019-2020.
