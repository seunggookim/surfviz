function [hc, contour] = tricontour(surf, val, levels, draw)
% plots linear interpolated contours on a 2D/3D triangulated suraface.
%
% [hc, contour]=tricontour(surf, val, levels)
%
% This function is based on FT_TRIPLOT() of the FIELDTRIP package
% (see http://www.ru.nl/neuroimaging/fieldtrip).
%
% Improvement is to form groups of edges so that MATLAB can handle them
% x180+ faster (e.g., 100k line obj [14 s] vs. 200 line obj, [0.075 s])
%
% Instead, now it takes about 3 sec to form line groups before drawing. 
% But this computation can be avoided by using the "contour" variable 
% itself. This is particularly useful if you need to draw same contours
% (i.e., isocurvature contours of the template surfaces) over and over 
% again for various views.
%
% This version is included in https://github.com/solleo/surfviz
%
% SEE ALSO: FSSS_VIEW, FT_TRIPLOT, FSSS_READ_ALL_FS_SURFS
%
% (cc) 2020, sgKIM, solleo@gmail.com

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
if ~exist('draw','var')
  draw = true;
end


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

tic
fprintf('Connecting edges...')

for cnt_indx = 1:length(levels)
  % ignore tiny differences from interpolation:
  p1 = round(contour(cnt_indx).intersect1*1000)/1000;
  p2 = round(contour(cnt_indx).intersect2*1000)/1000;
  verts = unique([p1;p2],'rows');
  [~, v1] = ismember(p1, verts, 'rows');
  [~, v2] = ismember(p2, verts, 'rows');
  edges = [v1, v2];
  
  isusededge = false(size(edges,1),1); % is this edge used or not?
  g = 1; % group index
  contour(cnt_indx).group = [];
  for i = 1:size(edges,1)
    % starting from the i-th edge: e_i = (v_i1, v_i2)
    if ~isusededge(i)
      isusededge(i) = true; % now the i-th edge is used.
      
      verts_thisgroup = edges(i,:);
      % put the vertices of the i-th edge to this group
      
      for v = 1:2 % consider v_i1 and v_i2 separately
        thisv = edges(i,v); % this vertex!
        
        newedges_thisv = find(~~sum(edges == thisv,2) & ~isusededge);
        % this is equivalent to:
        %    setdiff(Edges_including_this_vertex, Edges_used)
        % but much faster!
        
        while ~isempty(newedges_thisv)
          % looping until there's no further connection
          
          newv = edges(newedges_thisv(1),:);
          % just pick one edge from those new edges with "this vertex"
          
          newv(newv == thisv) = [];
          % find the new vertex
          
          if ~isempty(newv)
            % if you really have a new vertex
            
            if v == 1
              verts_thisgroup = [newv verts_thisgroup];
            else
              verts_thisgroup = [verts_thisgroup newv];
            end
            % THIS is the heart of this algorithm. We want to know whether
            % we should put this before or after the current list of
            % vertices to form a one curve.
            
            thisv = newv;
            % now let's look for further connections from the newly added
            % vertex by updating "thisv" with newv
            
            isusededge(newedges_thisv(1)) = true;
            % consider this edge used.
          else
            % this is possible if the edge is looping back (perhaps?)
            
            isusededge(newedges_thisv(1)) = true;
            % just consider this is used and move on.
          end
          
          newedges_thisv = find(~~sum(edges == thisv,2) & ~isusededge);
          % now update new edges with the updated isusededge
          
        end
      end
      contour(cnt_indx).group(g).xyz = verts(verts_thisgroup,:);
      g = g+1;
    end
  end
end
toc % 2.98 sec!

if draw
  tic
  fprintf('Drawing contours...');
  hc = hggroup;
  hold on
  for cnt_indx = 1:length(levels)
    for g = 1:numel(contour(cnt_indx).group)
      xyz = contour(cnt_indx).group(g).xyz;
      plot3(xyz(:,1),xyz(:,2),xyz(:,3), 'Parent',hc, 'linewidth',2);
    end
  end
  hold off
  toc % 0.049793
else
  hc = [];
end

end