# surfviz

MATLAB scripts for visualization of surface-mapped scalar values. Various FreeSurfer (http://freesurfer.net/) surface models (white, pial, inflated, sphere) can be used. OpenGL is required for transparent overlay. External MATLAB codes (`export_fig`, `brewermap`) are included. Largely inspired by visualization functions in `CAT12` (http://www.neuro.uni-jena.de/cat/).

## Install
Download/fetch the package and add all subdirectories to MATLAB path:
```
>>> addpath(genpath(DIRECTORY_WHERE_YOU_SAVED_FILES))
```

## Demo
```
subdir = fullfile(getenv('FREESURFER_HOME'),'subjects'); % Find FREESURFER default subject directory
surfs = fsss_read_all_FS_surfs('fsaverage',subdir); % Reading template surfaces
thns = cell(1,2);
thns{1} = read_curv(fullfile(subdir,'fsaverage','surf','lh.thickness')); % 'read_curv' is from FREESURFER MATLAB functions
thns{2} = read_curv(fullfile(subdir,'fsaverage','surf','rh.thickness'));
fsss_view(surfs, thns)
```
![](https://github.com/solleo/surfviz/blob/master/images/demo1.png)

You can also set the surface to visualize, layout, and a threshold via `cfg` structure:
```
cfg = struct('basesurf','PIAL','layout','1x2','thres',3);
fsss_view(surfs, thns, cfg)
```
![](https://github.com/solleo/surfviz/blob/master/images/demo2.png)

... and some captions and colorschemes: 
```
cfg = struct('basesurf','INFL','layout','2x2','thres',1,...
  'colorscheme','yellowblue',...
  'colorbartitle','Rel. Ctx. Thns.','colorbarxlabel','Z-socre');
fsss_view(surfs, thns_z, cfg)
```
![](https://github.com/solleo/surfviz/blob/master/images/demo4.png)

... finally you can create as many figures as you want without Matalb taking away your attention:
```
cfg = struct('basesurf','INFL','layout','2x2','thres',1,...
  'colorscheme','yellowblue',...
  'colorbartitle','Rel. Ctx. Thns.','colorbarxlabel','Z-socre',...
  'fname_png','demo4.png');
fsss_view(surfs, thns_z, cfg)

>>> ls demo4.png
demo4.png
```

See documentation for more information:
```
>>> doc fsss_view
```

(cc) Seung-Goo Kim
