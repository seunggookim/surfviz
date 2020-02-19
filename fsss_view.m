function [H, cfg] = fsss_view (surfs, data, cfg)
% FSSS_VIEW visualize data on surfs
%
% [USAGE]
% fsss_view (surfs, data)
% H = fsss_view (surfs, data)
% [H, cfg] = fsss_view (surfs, data, cfg)
%
%
% [INPUT]
% surfs {1x2} freesurfer surfaces structure read by FSSS_READ_ALL_FS_SURF
%
% data  {1x2} cells contain a vector (vertex-mapped scalar data) for each 
%             hemisphere. [NaN] will be displayed as transparent if openGL
%             is available.
%
% cfg   (1x1) structure to configure options (optional):
% -surfaces layout
%  .basesurf '1xN'  'INFL' (default) | 'WHITE' | 'PIAL'
%                   ('SMOOTHWM' | 'SPHERE' | 'SUPTEMP.*' if in surfs)
%  .layout   '1xN'  '2x4' (default) | '1x4' | '2x2' | '1x2' | '1x2oblq' |
%                   '1x2stmp'
%  .views    {1xN}  view angle, default is determined by layout:
%                   '2x4' {[-90 0],[90 0],[-90 0],[90 0],
%                          [-90 75],[90 -75],[-90 -75],[90 75]}
%                   '1x4' {[-90 0],[90 0],[-90 0],[90 0]}
%                   '2x2' {[-90 0],[90 0],[90 0],[-90 0]}
%                   '1x2' {[-90 0],[90 0]}
%                   '1x2oblq' {[-119 6],[106 -5]}
%                   '1x2stmp' {[0 90],[0 90]};
%  .hemis    [1xN]  1=left 2=right, default is determined by layout:
%                   '2x4' [1 1 2 2 1 1 2 2 ]
%                   '1x4' [1 1 2 2 ]
%                   '2x2' [1 2 1 2]
%                   '1x2' [1 2]
%  .figureposition [1x4] default is determined by layout:
%                   '2x4' [5   694   900   365]
%                   '1x4' [5   694   900   220]
%                   '2x2' [5   694   900   793]
%                   '1x2' [5   694   900   415]
%                   '1x2stmp' [5   694   900   615]
%  .figurecolor [1x3] | '1x1' default = 'k' if .thres ~= 0
%                                     = 'w' if .thres == 0
%
% -color schemes (presets)
%  .colorscheme    '1xN'  'darkspectral' | 'darkparula'
% -color coding (manual)
%  .basecolor '1xN'  'dark' (0.15/0.35) | 'bright' (0.7/1)
%  .facealpha [1x1]  default = 1 (overlay)
%  .colormap  [Nx3]  default = [bluish; gray; redish] if .thres ~= 0;
%                    brewermap(256,'spectral') if .thres == 0
%
% -- optional colormaps:
%   (.usetruncated) [1x1]  true | {false} for thresholded values
%  .caxis     [1x2] default is determined by signs of values:
%                   for positive values: [1%tile 99%tile]
%                   for sigend values: [-/+max(99%tile of abs)] 
%             alternatively, enter 'full' min/max of numerical values
%  .thres     [1x1] | [1x2]  default=0 (no thresholding)
%  .subthres  [1x1]  true | {false} - show subthreshold values with a=0.5
%  .masks     {1x2} cells contain vertex-mapped logical vectors for masking
%                   by default, non-cortical vertices are masked.
%
% -contour
%  .doisocurv       [1x1]  true | false (default) - surfs.WHITECURV
%  .clusid          {1x2}
%  .clusidcolor     [Cx3]  colormap for clusters, default = ones(C,3)
%  .clusidlinewidth [1x1]  linewidth, default = 0.5
%  .clusidlinestyle '1xN'  default = '-'
%
% -lighting
%  .camlight    [1x1]  {true} | false
%  .camposition [1x2]  default=[0 10]
%
% -histograms
%  .dohist        [1x1]  true (default) | false
%  .histfontsize  [1x1]  default=8
%  .histfontsize
%  .histxtick
%  .histxticklabel
%
% -colorbar
%  .docolorbar       [1x1]  {true} | false
%  .colorbartitle
%  .colorbarinterp   '1xN'  {'none'} | 'tex' | 'latex'
%  .colorbarfontsize
%  .colorbarxlabel
%  .colorbarxlabelinterp '1xN'  {'none'} | 'tex' | 'latex'
%
% -print (requires EXPORT_FIG)
%  .fname_png
%  .dpi
%
%
% [OUTPUT]
% H structure array (1xN) contains axes handles for:
% - surface axes
%  .axes
%  .basesurf
%  .oversurf
%  .light
%
% - histogram axes
%  .axes
%  .histfontsize
%  .histxtick
%  .histxticklabel
%
% - colorbar axes
%  .axes
%  .colorbar
%  .colorbarxtick
%  .colorbarxticklabel
%  .colorbartitle
%  .colorbarfontsize
%
% (cc) 2019, sgKIM. solleo@gmail.com
%
% See also VIEW_SURFDATA, FSSS_READ_ALL_FS_SURFS, EXPORT_FIG


%{
:TODO:
Style-schemes: "FSL" "FS" "unthres"

:HISTORY:
2019-12-17: removed "useparula" option from FSSS_VIEW and VIEW_SURFDATA
%}

%% C O N F I G ============================================================
if ~exist('cfg','var'), cfg=[]; end

%% -- Input dimension: forcing a column vector
for ihemi = 1:2
  if ~isvector(data{ihemi})
    error('data must contains a vector!')
  end
  data{ihemi} = reshape(data{ihemi},[],1);
end

%% -- Parallel workers: only Painter: can't handle tranparancy
if ~isempty(getCurrentWorker)
  warning('Parallel workers: only Painter: force unthresholded maps')
  cfg.basesurf = 'WHITE';
  cfg.nobasesurf = true;
%   if isfield(cfg,'thres'), cfg = rmfield(cfg, 'thres'); end
end

%% -- SUPTEMP
if isfield(cfg,'basesurf') && contains(cfg.basesurf,'SUPTEMP')
  if ~isfield(cfg,'layout')
    cfg.layout = '1x2stmp';
  end
  if ~isfield(cfg,'masks')
    cfg.masks = {surfs.ANNOT{1}.suptemp1_bin', surfs.ANNOT{2}.suptemp1_bin'};
  end
end

%% -- Color scheme
if isfield(cfg,'colorscheme')
  switch cfg.colorscheme
    case 'darkspectral'
      cfg.colormap = flipud(brewermap(256,'spectral'));
      cfg.figurecolor = 'k';
      cfg.basecolor = 'dark';
    case 'darkparula'
      cfg.colormap = parula(256);
      cfg.figurecolor = 'k';
      cfg.basecolor = 'dark';
    otherwise
      warning('%s is not valid colorscheme. Ignored', cfg.colorscheme)
  end
end
%% -- Layout
if ~isfield(cfg,'layout'), cfg.layout='2x4'; end
if ~isfield(cfg,'views')
  switch cfg.layout
    case '2x4'
      cfg.views = {...
        [-90 0],[90 0],[-90 0],[90 0], ...
        [-90 75],[90 -75],[-90 -75],[90 75]};
    case '1x4'
      cfg.views = {[-90 0],[90 0],[-90 0],[90 0]};
    case '2x2'
      cfg.views = {[-90 0],[90 0],[90 0],[-90 0]};
    case {'1x2big','1x2'}
      cfg.views = {[-90 0],[90 0]};
    case '1x2oblq'
      cfg.views = {[-120 10],[120 10]};
    case '1x2stmp'
      cfg.views = {[0 90],[0 90]};
  end
end
if ~isfield(cfg,'hemis')
  switch cfg.layout
    case '2x4'
      cfg.hemis = [1 1 2 2 1 1 2 2];
    case '1x4'
      cfg.hemis = [1 1 2 2];
    case '2x2'
      cfg.hemis = [1 2 1 2];
    case {'1x2','1x2big','1x2oblq','1x2stmp'}
      cfg.hemis = [1 2];
  end
end
if ~isfield(cfg,'axes')
  switch cfg.layout
    case '2x4'
      cfg.surfaxes = axeslayout([2 4],[.02 .02 0 0],[.02, .02, .0, .0]);
      cfg.surfaxes.y(1:4) = cfg.surfaxes.y(1:4)+0.02;
      cfg.surfaxes.y(5:8) = cfg.surfaxes.y(5:8)+0.15;
      cfg.histaxes = {[.05 .07 .25 .15],[.7 .07 .25 .15]};
      cfg.colorbaraxes = [.5-.15 .09 .3 .08];
    case '1x4'
      cfg.surfaxes = axeslayout([1 4],[.02 .02 0 0],[.02, .02, .0, .0]);
      cfg.surfaxes.y(1:4) = cfg.surfaxes.y(1:4)+0.15;
      cfg.histaxes = {[.05 .13 .25 .15],[.7 .13 .25 .15]};
      cfg.colorbaraxes = [.5-.15 .15 .3 .08];
    case '2x2'
      cfg.surfaxes = axeslayout([2 2],[.02 .02 0 0],[.02, .02, .0, .0]);
      cfg.surfaxes.y(1:2) = cfg.surfaxes.y(1:2)+0.05;
      cfg.surfaxes.y(3:4) = cfg.surfaxes.y(3:4)+0.15;
      cfg.histaxes = {[.05 .04 .25 .15],[.7 .04 .25 .15]};
      cfg.colorbaraxes = [.5-.15 .06 .3 .08];
    case '1x2'
      cfg.surfaxes = axeslayout([1 2],[.02 .02 0 0],[.02, .02, .0, .2]);
      cfg.surfaxes.y(1:2) = cfg.surfaxes.y(1:2)+0.06;
      cfg.histaxes = {[.08 .08 .23 .15],[.73 .08 .23 .15]};
      cfg.colorbaraxes = [.5-.15 .20 .3 .05];
    case {'1x2big','1x2oblq'}
      cfg.surfaxes = axeslayout([1 2],[.02 .02 0 0],[.02, .02, .0, .0]);
      cfg.surfaxes.y = cfg.surfaxes.y + 0.11;
      cfg.histaxes = {[.05 .07 .25 .15],[.7 .07 .25 .15]};
      cfg.colorbaraxes = [.5-.15 .09 .25 .05];
    case '1x2stmp'
      cfg.surfaxes = axeslayout([1 2],[.02 .02 0 0],[.02, .02, .0, .0]);
      cfg.surfaxes.y = cfg.surfaxes.y + 0.11;
      cfg.histaxes = {[.1 .05 .22 .10],[0.73 .05 .22 .10]};
      cfg.colorbaraxes = [.5-.25/2 .04 .25 .08];
  end
end
if ~isfield(cfg,'figureposition')
  switch cfg.layout
    case '2x4'
      cfg.figureposition = [5   694   900   365];
    case '1x4'
      cfg.figureposition = [5   694   900   220];
    case '1x2'
      cfg.figureposition = [5   694   450   220];
    case '2x2'
      cfg.figureposition = [5   694   900   793];
    case {'1x2big','1x2oblq'}
      cfg.figureposition = [5   694   900   415];
    case '1x2stmp'
      cfg.figureposition = [5   694   250   360];
  end
  if isunix
    cfg.figureposition(1) = 1922;
  end
end

%% -- Basesurf
if ~isfield(cfg,'basesurf'), cfg.basesurf = 'INFL'; end

%% -- Color range based on both hemisphere
if ~isfield(cfg,'masks')
  for hemi = 1:2
    cfg.masks{hemi} = true(size(data{hemi}));
  end
end
for hemi = 1:2
  cfg.masks{hemi} = cfg.masks{hemi} & surfs.ANNOT{hemi}.cortex;
end
if islogical(data{1,1}) || islogical(data{1,1})
  data = {uint8(data{1}), uint8(data{2})};
end
numvals = cat(1,data{:});
numvals(~cat(1,cfg.masks{:})) =[];
numvals(isnan(numvals)) = [];
numvals(isinf(numvals)) = [];
if isempty(numvals)
  warning('All vertices are masked')
  numvals = [0];
end
overallrange = [min(numvals) max(numvals)];

if ~isfield(cfg,'thres')
  cfg.thres = [0 0];
end
if numel(cfg.thres) == 1
  cfg.thres = [-abs(cfg.thres) abs(cfg.thres)];
end
if ~isfield(cfg,'caxis')
  if numel(unique(numvals)) < 2 % if it's binary values
    cfg.caxis = [0 1];
  elseif all(numvals>=0)
    cfg.caxis = [prctile(numvals,1) prctile(numvals,99)];
  else
    cfg.caxis = [-prctile(abs(numvals),99) +prctile(abs(numvals),99)];
  end
elseif ischar(cfg.caxis)
  switch cfg.caxis
    case 'full'
      cfg.caxis = [-max(abs(numvals)) max(abs(numvals))];
  end
end

%% -- Overlay values: integer or not?
if ~isfield(cfg,'isinterger')
  cfg.isinteger = isequal(numvals, round(numvals));
end

%% -- Initialize figure
figure
set(gcf, 'position', cfg.figureposition);
if isfield(cfg,'fname_png') % if fname_png is given, make it invisible
  set(gcf,'visible','off')
end


%% M A I N ================================================================
%% --- Plot surfaces
H = struct();
for iAxes = 1:numel(cfg.hemis)
  hemi = cfg.hemis(iAxes);
  axespos(cfg.surfaxes, iAxes);
  cfg.view = cfg.views{iAxes};
  % to show partial surfaces like supratemporal planes:
  idx = strfind(cfg.basesurf,'.');
  if ~isempty(idx)
    surf2show = ...
      surfs.(cfg.basesurf(1:idx-1)).(cfg.basesurf(idx+1:end)){hemi};
    data2show = data{hemi}(surf2show.vert_idx1_orig);
    cfg.mask = cfg.masks{hemi}(surf2show.vert_idx1_orig);
    cfg.curv = surfs.WHITECURV{hemi}(surf2show.vert_idx1_orig);
  else
    surf2show = surfs.(cfg.basesurf){hemi};
    data2show = data{hemi};
    cfg.mask = cfg.masks{hemi};
    cfg.curv = surfs.WHITECURV{hemi};
  end
  [h, cfg] = view_surfdata(surf2show, data2show, cfg);
  if cfg.isinteger
    h.oversurf.FaceColor = 'flat';
  end
  fn = fieldnames(h);
  for j=1:numel(fn), H(iAxes).(fn{j}) = h.(fn{j}); end
  %   axis on % for DEBUGGING
end

%% -- Plot histograms
if isequal(get(gcf,'color'), [1 1 1])
  cfg.histlinecolor = [.5 .5 .5];
  cfg.histfontcolor = [0 0 0];
elseif isequal(get(gcf,'color'), [0 0 0])
  cfg.histlinecolor = [.5 .5 .5];
  cfg.histfontcolor = [1 1 1];
end
if ~isfield(cfg,'histfontsize')
  if ismac
    cfg.histfontsize = 9;
  elseif isunix
    cfg.histfontsize = 8;
  elseif ispc
    cfg.histfontsize = 10;
  end
end

for hemi = 1:2
  iAxes = iAxes + 1;
  H(iAxes).axes = axes('position',cfg.histaxes{hemi});
  set(H(iAxes).axes,'color',get(gcf,'color'),...
    'xcolor',cfg.histfontcolor, 'ycolor',cfg.histfontcolor, ...
    'fontsize',cfg.histfontsize );
  hold on;
  numvals = cat(1,data{hemi});
  numvals(~cfg.masks{hemi}) = [];
  numvals(isnan(numvals)) = [];
  numvals(isinf(numvals)) = [];
  if numel(unique(numvals)) ~= 2
    numvals(numvals==0) = [];
  end
  if isempty(numvals)
    warning('[Hemi=%i] all vertices are masked. Skipping a histogram.',hemi);
    continue
  end
  
  numbins = min([200 numel(unique(numvals))]);
  if exist('histcounts','file')
    if cfg.isinteger
      [ci, edges] = histcounts(numvals, numbins, 'BinMethod','integers');
    else
      [ci, edges] = histcounts(numvals, numbins);
    end
    xi = 0.5 * (edges(1:end-1) + edges(2:end)); % center value of bins
  else
    [ci, xi] = hist(numvals, numbins);
  end
  
  % color-code bars:
  ind = zeroone(xi, cfg.caxis(1), cfg.caxis(2))*(size(cfg.colormap,1)-1);
  ind = 1+floor(ind);
  rgb = squeeze(ind2rgb(ind, cfg.colormap)) * 0.85;
  
  % suprathrs_neg:
  idx = xi <= cfg.thres(1);
  if sum(idx)
    h_neg = bar(xi(idx), ci(idx), 1, ...
      'edgecolor','none', 'facecolor','flat');
    h_neg.CData = rgb(idx,:);
  else
    h_neg = [];
  end
  
  % subthres:
  idx = cfg.thres(1) < xi & xi < cfg.thres(2);
  if sum(idx)
    h_subthres = bar(xi(idx), ci(idx), 1, ...
      'edgecolor','none', 'facecolor','flat');
    h_subthres.CData = rgb(idx,:);
  else
    h_subthres = [];
  end
  
  % suprathrs_pos:
  idx = cfg.thres(2) <= xi;
  if sum(idx)
    h_pos = bar(xi(idx), ci(idx), 1, ...
      'edgecolor','none', 'facecolor','flat');
    h_pos.CData = rgb(idx,:);
  else
    h_pos = [];
  end
  
  try
    xlim([edges(1) edges(end)])
  catch
    xlim(overallrange);
  end
  ylim0 = ylim;
  box on;
  
  % color-limits:
  line([cfg.caxis(1);cfg.caxis(1)],ylim0, ...
    'color', cfg.colormap(1,:)*.85, 'linewidth',0.5, 'linestyle',':');
  line([cfg.caxis(2);cfg.caxis(2)],ylim0, ...
    'color', cfg.colormap(end,:)*.85, 'linewidth',0.5, 'linestyle',':');
  
  % thresholding:
  if ~isequal(cfg.thres, [0 0])
    line([cfg.thres;cfg.thres],[ylim0' ylim0'],...
      'color',cfg.histlinecolor, 'linewidth',0.5, 'linestyle',':');
    idx = cfg.thres(1) < xi & xi < cfg.thres(2);
    if any(idx)
      if ~cfg.subthres
        h_subthres.CData = repmat(cfg.basecolormap(2,:),[sum(idx),1]);
      else
        h_subthres.FaceAlpha = 0.5;
      end
    end
  end
  
  if isfield(cfg, 'histxtick')
    set(H(iAxes).axes, 'xtick', cfg.histxtick)
  end
  if isfield(cfg, 'histxticklabel')
    set(H(iAxes).axes, 'xticklabel', cfg.histxticklabel)
  end
  H(iAxes).hist = [h_neg, h_subthres, h_pos];
end

%% -- Plot colorbar
numvals = cat(1,data{:});
numvals(~cat(1,cfg.masks{:})) =[];
numvals(isnan(numvals)) = [];
numvals(isinf(numvals)) = [];
if numel(unique(numvals)) ~= 2
  numvals(numvals==0) = [];
end
if ~isempty(numvals) %&& ~strcmp(cfg.layout,'1x2')
  iAxes = iAxes + 1;
  H(iAxes).axes = axes('position',cfg.colorbaraxes);
  if ~isfield(cfg,'colorbarfontsize')
    if ismac
      cfg.colorbarfontsize = 13;
    elseif isunix
      cfg.colorbarfontsize = 10;
    elseif ispc
      cfg.colorbarfontsize = 11;
    end
  end
  H(iAxes).colorbar = colorbar(H(iAxes).axes, 'location','north', ...
    'fontsize',cfg.colorbarfontsize*0.9);
  caxis(H(iAxes).axes, cfg.caxis);
  if cfg.isinteger && (numel(numvals)<50) % descritize colorbar for integer values
    numbins = min([200 numel(unique(numvals))]);
    if exist('histcounts','file')
      [~, edges] = histcounts(numvals, numbins, 'BinMethod','integers');
      xi = 0.5 * (edges(1:end-1) + edges(2:end)); % center value of bins
    else
      [~, xi] = hist(numvals, numbins);
    end
    ind = zeroone(xi, cfg.caxis(1), cfg.caxis(2))*(size(cfg.colormap,1)-1);
    ind = 1+floor(ind);
    cfg.colormap = squeeze(ind2rgb(ind, cfg.colormap));
  end
  colormap(H(iAxes).axes, cfg.colormap)
  axis off;
  if ~isfield(cfg,'colorbarxlabelinterp')
    cfg.colorbarxlabelinterp = 'none';
  end
  if isfield(cfg,'colorbarxlabel')
    xlabel(H(iAxes).colorbar, cfg.colorbarxlabel, ...
      'fontsize', cfg.colorbarfontsize, 'color',cfg.histfontcolor, ...
      'interpreter',cfg.colorbarxlabelinterp);
  end
  if ~isfield(cfg,'colorbarinterp')
    cfg.colorbarinterp='none';
  end
  if isfield(cfg,'colorbartitle')
    title(H(iAxes).colorbar, cfg.colorbartitle, ...
      'fontsize', cfg.colorbarfontsize*1.2,'fontweight','bold',...
      'color',cfg.histfontcolor, 'interp',cfg.colorbarinterp);
  end
  if ~isfield(cfg,'colorbarxtick')
    xtick = get(H(iAxes).colorbar, 'xtick');
    if (cfg.isinteger) && (numel(numvals)<50)
      warning('Would this work?')
      step = (range(xtick)-1)/range(xtick);
      xtick = [xtick(1)+step:step:xtick(2)-step];
    end
    cfg.colorbarxtick = xtick;
  end
  set(H(iAxes).colorbar,'color',get(gcf,'color'),...
    'xcolor',cfg.histfontcolor, 'ycolor',cfg.histfontcolor,...
    'xtick',cfg.colorbarxtick);
  if isfield(cfg,'colorbarxticklabel')
    set(H(iAxes).colorbar, 'xticklabel', cfg.colorbarxticklabel)
  end
end

%% -- Print
% set(gcf,'visible','on')
if isfield(cfg,'fname_png')
  if ~isfield(cfg,'dpi')
    cfg.dpi = 300;
  end
  if ~isempty(getCurrentWorker)
    rendopt = '-painters';
  else
    rendopt = '-opengl';
  end
  export_fig(cfg.fname_png,['-r',num2str(cfg.dpi)],rendopt)
  close(gcf)
end

if ~nargout, clear H cfg; end

end




