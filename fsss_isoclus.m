function surfs = fsss_isoclus(surfs, label, cfg)
% surfs = fsss_isoclus(surfs, label, cfg)
%
% taking `label` cell array <1x2> as input
% returning surfs.(cfg.basesurf).isoclus structure to be used in FSSS_VIEW
%
% (cc) 2020, sgKIM.

if ~exist('cfg','var')
  cfg = [];
end
if ~isfield(cfg,'basesurf')
  cfg.basesurf = 'INFL';
end
for s = 1:2
  [~,~,idx] = unique(label{s}); % convert it into continous labels from 1 to N
  surfs.(cfg.basesurf){s}.isoclus = struct();
  for i = 1:max(idx)
    smoothclus = SurfStatSmooth(...
      double(idx==i)', surfs.(cfg.basesurf){s}, 4); 
    % smoothing a little to beutify it a little (but now you need SurfStat)
    
    [~,c] = tricontour(surfs.(cfg.basesurf){s}, smoothclus, 0.5, false);
    surfs.(cfg.basesurf){s}.isoclus(i).group = c.group;
  end
end