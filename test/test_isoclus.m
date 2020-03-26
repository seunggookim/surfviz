clear
[v,f] = FS_read_surf('/Applications/freesurfer/subjects/fsaverage/surf/lh.pial_semi_inflated');
lh.vertices = v;
lh.faces = f;
curv = read_curv('/Applications/freesurfer/subjects/fsaverage/surf/lh.curv');
[~,c] = tricontour(lh, curv, 0, false);
lh.isocurv = c.group;
thns = read_curv('/Applications/freesurfer/subjects/fsaverage/surf/lh.thickness');
%%
view_surfdata(lh, thns, struct('basecolormap',[0 0 0; 1 1 1], 'colormap',jet, 'doisocurv',1,'figurecolor','w'))

%%
fname = '/Users/sol/fsaverage/label/lh.HCP-MMP1.annot';
[verts, label, cot] = read_annotation( fname, 0 );
%%
clf
view_surfdata(lh, label, struct('basecolormap',[0 0 0; 1 1 1], 'colormap',cot.table(:,1:3)/256, 'doisocurv',1,'figurecolor','w'))

%%
clear
[v,f] = FS_read_surf('/Applications/freesurfer/subjects/fsaverage/surf/lh.inflated');
lh.vertices = v;
lh.faces = f;
curv = read_curv('/Applications/freesurfer/subjects/fsaverage/surf/lh.curv');
[~,c] = tricontour(lh, curv, 0, false);
lh.isocurv = c.group;
thns = read_curv('/Applications/freesurfer/subjects/fsaverage/surf/lh.thickness');
fname = '/Users/sol/fsaverage/label/lh.HCP-MMP1.annot';
[verts, label, cot] = read_annotation( fname, 0 );
clf
view_surfdata(lh, label, struct('basecolormap',[0 0 0; 1 1 1], 'colormap',cot.table(:,1:3)/256, 'doisocurv',1,'figurecolor','w'))
h = findobj('type','patch');
h(1).FaceColor = 'flat';
%%
smoothlh = lh;
%%
[~,~,idx] = unique(label);
lh.isoclus = struct();
for i = 1:max(idx)
  smoothclus = SurfStatSmooth(double(idx==i)', lh, 4);
  [~,c] = tricontour(lh, smoothclus, 0.5, false);
  lh.isoclus(i).group = c.group;
end
%%
clf
view_surfdata(lh, thns, struct('basecolormap',[0 0 0; 1 1 1], 'caxis',[1 4], 'colormap',flipud(brewermap(256,'spectral')), 'doisocurv',1,'figurecolor','w','doisoclus',1,'isocluslinewidth',2,'isocluscolor',cot.table(:,1:3)/256))
%%
figure
view_surfdata(lh, thns, struct('basecolormap',[0 0 0; 1 1 1], 'caxis',[1 4], 'colormap',gray, 'doisocurv',0,'figurecolor','w','doisoclus',1,'isocluslinewidth',2,'isocluscolor',brewermap(max(idx),'Dark2')))
%%
figure
view_surfdata(lh, thns, struct('basecolormap',[0 0 0; 1 1 1], 'caxis',[1 4], 'colormap',gray, 'doisocurv',0,'figurecolor','w','doisoclus',1,'isocluslinewidth',2,'isocluscolor',cot.table(:,1:3)/256))
%%
fsavg = fsss_read_all_FS_surfs('fsaverage');
labels = {};
[~, labels{1}, cot] = read_annotation( '/Users/sol/fsaverage/label/lh.HCP-MMP1.annot', 0 );
[~, labels{2}, cot] = read_annotation( '/Users/sol/fsaverage/label/rh.HCP-MMP1.annot', 0 );
fsavg = fsss_isoclus(fsavg, labels); % this takes a while for granular annotations
%%
gray256=gray(256);
fsss_view(fsavg, fsavg.THNS, struct('layout','2x4','colormap',gray256(65:end,:), 'isocluslinewidth',1,'isocluscolor',cot.table(:,1:3)/256))
%%

