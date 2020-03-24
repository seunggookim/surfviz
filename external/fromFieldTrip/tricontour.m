function [h]=tricontour(surf, val, levels)
% plots linear interpolated contours on a 2D/3D triangulated suraface.
%
% This function is from FT_TRIPLOT() of the FIELDTRIP package
% (see http://www.ru.nl/neuroimaging/fieldtrip).
%
% tricontour(pnt, tri, val)

% Copyright (C) 2001=2006, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: triplot.m 6321 2012-08-03 16:02:37Z roevdmei $

pnt = surf.vertices;
tri = surf.faces;
absmax = max(abs([min(val) max(val)]));
if ~exist('levels','var')
  levels = linspace(-absmax,absmax,21);
end
% levels = 1e-10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute contours for 2D or 3D triangulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Computing contour..')
tic
% map values onto vertex list
triangle_val = val(tri);
% find min value for each triangle
triangle_min = min(triangle_val, [], 2);
% find max value for each triangle
triangle_max = max(triangle_val, [], 2);

for cnt_indx=1:length(levels)
  cnt = levels(cnt_indx);
  % find triangles with [min <= current level <= max]
  use = cnt>=triangle_min & cnt<=triangle_max;
  % initialize
  counter = 0;
  intersect1 = [];
  intersect2 = [];
  
  for tri_indx=find(use)' % find the triangles
    pos  = pnt(tri(tri_indx,:), :);
    v(1) = triangle_val(tri_indx,1);
    v(2) = triangle_val(tri_indx,2);
    v(3) = triangle_val(tri_indx,3);
    la(1) = (cnt-v(1)) / (v(2)-v(1)); % abcissa between vertex 1 and 2
    la(2) = (cnt-v(2)) / (v(3)-v(2)); % abcissa between vertex 2 and 3
    la(3) = (cnt-v(3)) / (v(1)-v(3)); % abcissa between vertex 1 and 2
    abc(1,:) = pos(1,:) + la(1) * (pos(2,:) - pos(1,:));
    abc(2,:) = pos(2,:) + la(2) * (pos(3,:) - pos(2,:));
    abc(3,:) = pos(3,:) + la(3) * (pos(1,:) - pos(3,:));
    counter = counter + 1;
    sel     = find(la>=0 & la<=1);
    intersect1(counter, :) = abc(sel(1),:);
    intersect2(counter, :) = abc(sel(2),:);
  end
  
  % remember the details for external reference
  contour(cnt_indx).level = cnt;
  contour(cnt_indx).n     = counter;
  contour(cnt_indx).intersect1 = intersect1;
  contour(cnt_indx).intersect2 = intersect2;
end

% collect all different contourlevels for plotting
intersect1 = [];
intersect2 = [];
cntlevel   = [];
for cnt_indx=1:length(levels)
  intersect1 = [intersect1; contour(cnt_indx).intersect1];
  intersect2 = [intersect2; contour(cnt_indx).intersect2];
  cntlevel   = [cntlevel; ones(contour(cnt_indx).n,1) * levels(cnt_indx)];
end

X = [intersect1(:,1) intersect2(:,1)]';
Y = [intersect1(:,2) intersect2(:,2)]';
C = [cntlevel(:)     cntlevel(:)]';
if size(pnt,2)>2
  Z = [intersect1(:,3) intersect2(:,3)]';
else
  Z = zeros(2, length(cntlevel));
end
toc

% tic
% fprintf('Connecting edges...')
% for cnt_indx = 1:length(levels)
%   p1 = (contour(cnt_indx).intersect2);
%   p2 = (contour(cnt_indx).intersect1);
%   igrp = 1;
%   G = cell(1,size(p1,1));
%   G{igrp} = [p1(1,:)];
%   for j = 2:size(p1,1)
%     G{igrp} = [G{igrp}; p1(j,:); p2(j,:)];
%     if any(p1(j,:) ~= p2(j-1,:))
%       igrp = igrp + 1;
%     end
%   end
%   G(igrp:end) = [];
%   igrp
%   for igrp = 1:numel(G)
%     h = plot3(G{igrp}(:,1), G{igrp}(:,2), G{igrp}(:,3),'color','k');
%   end
% %   tic
% %   verts = unique([p1;p2],'rows');
% %   [~,v1] = ismember(p1, verts, 'rows');
% %   [~,v2] = ismember(p2, verts, 'rows');
% end
% toc % 0.04 sec?

tic
fprintf('Connecting edges...')
for cnt_indx = 1:length(levels)
  p1 = (contour(cnt_indx).intersect1);
  p2 = (contour(cnt_indx).intersect2);
  tic
  verts = unique([p1;p2],'rows');
  [~,v1] = ismember(p1, verts, 'rows');
  [~,v2] = ismember(p2, verts, 'rows');
  edges = [v1 v2];

  edgegroup = nan(size(edges,1),1);
  %   group_index = 1;
  for i = 1:size(edges,1)
    for v = 1:2
      thisv = edges(i,v);
      edges_including_thisv = ~~sum(edges == thisv,2);
      if all(isnan(edgegroup(edges_including_thisv)))
        edgegroup(edges_including_thisv) = i;
      else
        edgegroup(edges_including_thisv) = nanmin(edgegroup(edges_including_thisv));
      end
    end
  end

  group_indices = unique(edgegroup(~isnan(edgegroup)))';
  for g = group_indices
    thisedges = unique(edges(edgegroup == g,:),'rows');
    


    plot3(verts(edgegroup == g,1), verts(edgegroup == g,2), verts(edgegroup == g,3), ...
      'color','k')
  end

  toc

end
toc % 0.04 sec?

% disp('rendering')
% tic;
% % make black-white contours
% hc = zeros(size(cntlevel));
% for i=1:length(cntlevel)
%   if cntlevel(i)>0
%     linestyle = '-';
%     linewidth = 1;
%   elseif cntlevel(i)<0
%     linestyle = '--';
%     linewidth = 1;
%   else
%     linestyle = '-';
%     linewidth = 2;
%   end
%   h1 = patch('XData', X(:,i), 'Ydata', Y(:,i), ...
%     'ZData', Z(:,i), 'CData', C(:,i)*nan, ...
%     'facecolor','none','edgecolor','black', ...
%     'linestyle', linestyle, 'linewidth', linewidth, ...
%     'userdata',cntlevel(i));
%   hc(i) = h1;
% end
% toc % 4 sec?!


% tic
% fprintf('Drawing contours...');
% h = line(X,Y,Z, 'color','k','linewidth',2);
% toc % PROBLEM IS THAT NOW MATLAB CREATES TOO MANY OBJECTS SINCE 2006!

% tic
% fprintf('Drawing contours...');
% h = plot3(X(:),Y(:),Z(:), 'color','k');
% toc

disp([])

end