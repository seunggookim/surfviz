clear
surfs = fsss_read_all_FS_surfs('fsaverage');
%%
clf
hc = tricontour(surfs.INFL{1}, surfs.WHITECURV{1}, 1e-10)
axis image
% 2.98 sec
% %%
% ft_triplot(surfs.INFL{1}.vertices, surfs.INFL{1}.faces, surfs.WHITECURV{1},...
%   'contour_BW',1e-10)
% % 8.74 sec
%%
clear
surfs = fsss_read_all_FS_surfs('fsaverage6',[], struct('isocurv','INFL'));
%%
fsss_view(surfs, surfs.THNS, struct('thres',2, 'isocurvcolor',[1 1 1]))
%%
fsss_view(surfs, surfs.THNS, struct('isocurvcolor',[1 1 1]))
