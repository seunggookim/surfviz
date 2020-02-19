# surfviz

MATLAB scripts for visualization of surface-mapped scalar values. Various FreeSurfer (http://freesurfer.net/) surface models (white, pial, inflated, sphere) can be used. OpenGL is required for transparent overlay. External MATLAB codes (`export_fig`, `brewermap`) are included.

## Install
Download/fetch the package and add all subdirectories to MATLAB path:
```
>>> addpath(genpath(DIRECTORY_WHERE_YOU_SAVED_FILES))
```

## Demo
```
surfs = fsss_read_all_FS_surfs('bert',getenv('SUBJECTS_DIR')); % Reading surfaces of FREESURFER's example subject "bert" assuming the evironmental variable $SUBJECT_DIR is set as $FREESURFER_HOME/subjects
thns = cell(1,2);
thns{1} = read_curv(fullfile(getenv('SUBJECTS_DIR'),'bert','surf','lh.thickness')); % 'read_curv' is from FREESURFER MATLAB functions
thns{2} = read_curv(fullfile(getenv('SUBJECTS_DIR'),'bert','surf','rh.thickness'));
fsss_view(surfs, thns)
```
(cc) Seung-Goo Kim
