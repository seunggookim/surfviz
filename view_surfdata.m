function [H, cfg] = view_surfdata (surf, data, cfg)
% [H, cfg] = VIEW_SURFDATA (surf, data, cfg)
%
% visualize scaled data on a surface
%
% [INPUT]
%  .surf     (1x1)  MATLAB patch structure
%  .data     [Vx1]  Vertex-mapped scalar data
%  .view     [1x2]  view angle, default = [-90 0]
%
% -color schemes (presets)
%  .style    '1xN'  'spectral' | 'spectraltrans' | 'fsl' | 'freesurfer' 
% -color coding (manual)
%  .basecolor '1xN'  'dark' (0.15/0.35) | 'bright' (0.7/1)
%  .colormap  [Nx3]  default = [bluish; gray; redish] if .thres ~= 0; 
%                    brewermap(256,'spectral') if .thres == 0
% -- optional colormaps:
%  .caxis     [1x2]  default = [1%tile 99%tile] for positive values; 
%                    [-/+max(99%tile of abs)] for sigend values
%  .thres     [1x1] | [1x2]  default=0
%  .subthres  [1x1]  true | false (default) - show subthreshold values
%                    transparently
%  .mask      [Vx1]  vertex-mapped logical vectors for masking
%
% -contours
%  .doisocurv       [1x1]  true (surfs.WHITECURV) | false (default) 
%  .clusid          {1x2}  
%  .clusidcolor     [Cx3]  colormap for clusters, default = ones(C,3)
%  .clusidlinewidth [1x1]  linewidth, default = 0.5
%  .clusidlinestyle '1xN'  default = '-'
%
% -lighting
%  .camlight    [1x1]  true (default) | false
%  .camposition [1x2]  default = [0 10]
%
%  .figurecolor [1x3] | '1x1' default = 'k' if .thres ~= 0
%                                     = 'w' if .thres == 0
%
% [OUTPUT]
% H.basesurf
% H.oversurf
% H.light
% H.axes
%
%
% See also VIEWSURFDATA
%
% (cc) 2019, sgKIM. mailto://solleo@gmail.com   https://ggooo.wordpress.com

%--------------------------------------------------------------------------
if ~exist('cfg','var'), cfg=[]; end
% for SurfStat surfaces
if isfield(surf,'tri') && ~isfield(surf,'faces')
 surf.faces = single(surf.tri);
end
if isfield(surf,'coord') && ~isfield(surf,'vertices')
 surf.vertices = single(surf.coord)';
end

%% Color axis                                                              
if islogical(data), data = uint8(data); end
numvals = data;
if isfield(cfg,'mask'), numvals(~cfg.mask) =[]; end
numvals(isnan(numvals)) = [];
numvals(isinf(numvals)) = [];

if ~isfield(cfg,'thres')
  cfg.thres = [0 0];
end
if numel(cfg.thres) == 1
  cfg.thres = [-abs(cfg.thres) abs(cfg.thres)];
end
if ~isfield(cfg,'caxis')
  if numel(unique(numvals)) < 2 % if it's binary values
    cfg.caxis = [min(numvals) max(numvals)];
  else % Winsorization
    if all(numvals>=0)
      cfg.caxis = [prctile(numvals,1) prctile(numvals,99)];
    else
      cfg.caxis = [-prctile(abs(numvals),99) +prctile(abs(numvals),99)];
    end
  end
end

% %% Color schemes
% NUMCOLORLEVELS = 256;
% if ~isfield(cfg,'subthres')
%   cfg.subthres = 0;
% end
% if ~isfield(cfg,'colormap')
%   if isequal(cfg.thres, [0 0])
%     cfg.colormap = flipud(brewermap(NUMCOLORLEVELS,'spectral'));
%     cfg.figurecolor = 'w';
%     cfg.basecolor = 'bright';
%   else % thresholded
%     ind = linspace(cfg.caxis(1),cfg.caxis(2),NUMCOLORLEVELS);
%     NUMLEVELNEG = sum(ind <= cfg.thres(1));
%     NUMLEVELSUBTHR = sum(cfg.thres(1) < ind & ind < cfg.thres(2));
%     NUMLEVELPOS = sum(cfg.thres(2) <= ind);
%     cfg.colormap = [
%       flipud(sgcolormap('DBC3',NUMLEVELNEG));
%       0.35*ones(NUMLEVELSUBTHR,3);
%       sgcolormap('DRY',NUMLEVELPOS)];
%     cfg.figurecolor = 'k';
%     cfg.basecolor = 'dark';
%   end
% else
%   if ~isfield(cfg,'figurecolor')
%     cfg.figurecolor = 'w';
%   end
%   if ~isfield(cfg,'basecolor')
%     cfg.basecolor = 'bright';
%   end
% end
% switch cfg.basecolor
%   case 'dark'
%     cfg.basecolormap = [.15 .15 .15; .35 .35 .35];
%   case 'bright'
%     cfg.basecolormap = [.6 .6 .6; .8 .8 .8];
% end
% if cfg.subthres
%   cfg.basecolormap = [.7 .7 .7; .9 .9 .9];
%   cfg.colormap = flipud(brewermap(NUMCOLORLEVELS,'spectral'));
%   cfg.facealpha = 0.6;
%   cfg.figurecolor = 'k';
% end

%% base image: binary curvature ([0,1]) or just [1])
if ~isfield(cfg,'curv')
  img_base = ones(size(data));
else
  img_base = single(cfg.curv<0);
end
if cfg.subthres
  img_base = single(cfg.curv<0);
  img_base(cfg.mask) = ones(sum(cfg.mask),1);
end
img_base_rgb = squeeze(ind2rgb(img_base+1, cfg.basecolormap)); 
% see CAT_SURF_RENDER for an example of using ind2rgb

%% overlay image
img_over = data;
img_over(cfg.thres(1)<data & data<cfg.thres(2)) = nan; % thresholding
if isfield(cfg,'mask'), img_over(~cfg.mask) = nan; end % masking
if ~isempty(getCurrentWorker)
  img_over(isnan(img_over)) = 0;
end

%% Surfaces
H.axes = gca;
hold on
V = surf.vertices;
F = surf.faces;
if isfield(cfg,'nobasesurf') && cfg.nobasesurf
  H.basesurf = [];
else
  H.basesurf = patch('faces',F, 'vertices',V, 'facecolor','interp', ...
    'edgecolor', 'none', 'FaceVertexCData', img_base_rgb, ...
    'ambientstrength',0.4, 'diffusestrength',0.8, 'specularstrength',0);
end
% see CAT_SURF_RENDER for a good exmaple of PATCH options
if ~isfield(cfg,'facealpha'), cfg.facealpha = 1; end
H.oversurf = patch('faces',F, 'vertices',V, 'facecolor','interp', ....
  'edgecolor', 'none', 'FaceVertexCData', img_over, ...
  'facealpha', cfg.facealpha, ...
  'ambientstrength',0.4, 'diffusestrength',0.8, 'specularstrength',0);

%% subthreshold overlay image
if cfg.subthres
  img_over = data;
  if isfield(cfg,'mask'), img_over(~cfg.mask) = nan; end % masking
  H.subthsurf = patch('faces',F, 'vertices',V, 'facecolor','interp', ...
  'edgecolor', 'none', 'FaceVertexCData', img_over, 'facealpha',0.5, ...
  'ambientstrength',0.4, 'diffusestrength',0.8, 'specularstrength',0);
end
axis off, axis tight, axis image
hold off
if isfield(cfg, 'facecolor')
  set(H.oversurf, 'FaceColor',cfg.facecolor); 
end
if ~isfield(cfg, 'view'),  cfg.view = [-90 0]; end
view(cfg.view)
set(gca, 'colormap',cfg.colormap)
caxis(gca, cfg.caxis);
set(gcf,'color',cfg.figurecolor)

%% Lighting
if ~isfield(cfg,'camlight'), cfg.camlight = 1; end
if ~isfield(cfg,'camposition'), cfg.camposition=[0 10]; end  
if cfg.camlight
 if ismac && isempty(getCurrentWorker) % MAC & not by a worker
   %===lines from CAT_SURF_RENDER (CAT12/SPM toolbox)
   % REF: 
   % set inner light 
   H.light = light('Position',[0 0 0]);
   set(H.basesurf,'BackFaceLighting','unlit');
   set(H.oversurf,'BackFaceLighting','unlit');
   if isfield(H,'subthsurf')
     set(H.subthsurf,'BackFaceLighting','unlit');
   end
   %===lines from CAT_SURF_RENDER
 else
   H.light = camlight(cfg.camposition(1), cfg.camposition(2));
 end
end

%% Cluster contours
if isfield(cfg,'clusid')
  if ~isfield(cfg,'clusidlinewidth')
    cfg.clusidlinewidth = 0.5;
  end
  if ~isfield(cfg,'clusidlinestyle')
    cfg.clusidlinestyle = '-';
  end
  if ~isfield(cfg,'clusidcolor')
    cfg.clusidcolor = ones(3,max(cfg.clusid));
  end
  for c=1:max(cfg.clusid) % for each cluster ID
    [~,H.contour_clusid(c)] = ft_triplot(V, F, double(cfg.clusid==c), ...
      'contour_bw', 1e-10);
    set(H.contour_clusid(c), 'color',cfg.clusidcolor(c,:), ...
      'linewidth',cfg.clusidlinewidth, 'linestyle',cfg.clusidlinestyle);
  end
end

%% Isocurvature contours
if ~isfield(cfg,'doisocurv')
  cfg.doisocurv = 0;
end
if cfg.doisocurv
 if ~isfield(cfg,'isocurvcolor')
  cfg.isocurvcolor=[0 0 0];
 end
 if ~isfield(cfg,'isocurvcontourwidth')
  cfg.isocurvcontourwidth=0.5;
 end
 if ~isfield(cfg,'noisocurvcontour')
  hold on;
  [~,H.contour_isocurv] = ft_triplot(V, F, cfg.curv, 'contour_bw', 1e-10);
  set(H.contour_isocurv, 'color',cfg.isocurvcolor, ...
    'linewidth', cfg.isocurvcontourwidth);
 end
end
end
