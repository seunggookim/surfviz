function fsss_view1 (y, side, subj, fsdir, doIsoCurv)
% fsss_view1 (y, side, subj, fsdir, doIsoCurv)

if nargin==0
help fsss_view1
return
end
%hemi={'lh','rh'};
if ~exist('subj','var') || isempty(subj)
subj='fsaverage6'; 
end
if ~exist('side','var') || isempty(side)
side=1; 
end
if ~exist('fsdir','var') || isempty(fsdir)
fsdir='/scr/vatikan1/bin/freesurfer-Linux-centos6_x86_64-stable-pub-v5.3.0//subjects/';
end

surfs = fsss_read_all_FS_surfs(subj,fsdir);
if side == 2
V=[90, 0];
else
V=[-90 0];
end
cfg=struct('view',V);
if ~exist('doIsoCurv','var') || isempty(doIsoCurv)
else
cfg.curv = surfs.WHITECURV{side};
end

if sum(~~y(:)) < size(surfs.INFL{side}.vertices,1)*0.5;
val=y;
y = double(surfs.WHITECURV{side} > 0);
y(~~val) = val(~~val)+1;
num_level=numel(unique(val(~~val)));

cfg.colormap=[.8 .8 .8; .99 .99 .99; hsv]; 
%cfg.colormap=[.8 .8 .8; .95 .95 .95; 1 0 0; 0 1 0; 0 0 1; 0 1 1]; 
%cfg.caxis=[0 num_level+1];
end

view_trisurf(surfs.INFL{side}, y, cfg)
end
