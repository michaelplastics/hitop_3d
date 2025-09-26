function viewSkeletonVsVoxel(segDir, voxMatPath)
% viewSkeletonVsVoxel(segDir, voxMatPath)
% 1) Aligns skeleton segments (OBJ polylines) to voxels.
% 2) Maps each solid voxel to the nearest skeleton segment.
% 3) Opens two displays:
%    • Alignment (as before)
%    • Solid, cubical voxel clusters (per segment) — INTERACTIVE:
%        - Hover a cluster to highlight
%        - Double-click a cluster to assign an integer label
%        - Save (button or 'S') selections to JSON (saved one dir ABOVE segDir)
%        - Clear ('C') selections
%
% JSON includes: segment_index, voxel_count, label, and voxel coordinates [x y z].

% ---------- Load voxel volume ----------
S = load(voxMatPath);
voxBin = get_voxbin_from_struct(S);         % [ny, nz, nx] logical
[ny,nz,nx] = size(voxBin);

% Solid voxels (index space, x-y-z)
[idxY, idxZ, idxX] = ind2sub(size(voxBin), find(voxBin));
voxCoord = [idxX idxY idxZ];

% ---------- Load skeleton segments ----------
objs = dir(fullfile(segDir, '*.obj'));
assert(~isempty(objs), 'No OBJ segments found in: %s', segDir);

Vseg = cell(numel(objs),1);
for i = 1:numel(objs)
    Vseg{i} = load_obj_vertices(fullfile(objs(i).folder, objs(i).name));
end
nSeg = numel(Vseg);

% ---------- Auto-align skeleton to voxel frame ----------
[R, t] = choose_alignment(voxBin, Vseg);
for i = 1:nSeg
    Vseg{i} = (R * Vseg{i}.' ).' + t;
end

% ========== Display 1: Alignment ==========
fig1 = figure('Color','w','Name','Alignment: Skeleton vs Voxel');
ax1  = axes(fig1); hold(ax1,'on'); axis(ax1,'equal','off'); set(ax1,'YDir','normal'); view(ax1,3);
try, enableDefaultInteractivity(ax1); end

voxXYZ = [idxX idxY idxZ];
step = max(1, floor(numel(idxX)/30000));
scatter3(ax1, voxXYZ(1:step:end,1), voxXYZ(1:step:end,2), voxXYZ(1:step:end,3), 4, [0.85 0.85 0.85], '.');

for i = 1:nSeg
    P = Vseg{i};
    if size(P,1) >= 2
        plot3(ax1, P(:,1), P(:,2), P(:,3), 'LineWidth', 1.5);
    else
        plot3(ax1, P(:,1), P(:,2), P(:,3), 'o', 'MarkerSize', 4);
    end
end
title(ax1, sprintf('Alignment  (%dx%dx%d)', nx, ny, nz));

% ========== Mapping: voxel → nearest segment ==========
[segVol, vox2seg] = map_voxels_to_segments(voxBin, voxCoord, Vseg);

% Precompute per-segment voxel coordinate lists [x y z] for JSON
segVoxelCoords = cell(1, nSeg);
for s = 1:nSeg
    if any(segVol{s}(:))
        [yy,zz,xx] = ind2sub(size(voxBin), find(segVol{s}));
        segVoxelCoords{s} = [xx(:) yy(:) zz(:)];
    else
        segVoxelCoords{s} = zeros(0,3);
    end
end

% ========== Display 2: Solid, cubical, no-shadows (interactive) ==========
fig2 = figure('Color','w','Name','Voxel→Segment Mapping (cubical, interactive)');
ax2  = axes(fig2);
hold(ax2,'on'); axis(ax2,'equal','off'); set(ax2,'YDir','normal'); view(ax2,3);
box(ax2,'on'); daspect(ax2,[1 1 1]);
lighting(ax2,'none'); shading(ax2,'flat'); material(ax2,'dull');
set(fig2, 'Renderer','opengl');
try, enableDefaultInteractivity(ax2); end  % pan/zoom/rotate enabled

C = simple_distinguishable_colors(nSeg);
Hseg = gobjects(nSeg,1);
for s = 1:nSeg
    mask = segVol{s};
    if ~any(mask(:)), continue; end
    [Vquad, Fquad] = build_voxel_mesh_quads(mask);
    Hseg(s) = patch(ax2, 'Vertices', Vquad, 'Faces', Fquad, ...
        'FaceColor', C(s,:), ...
        'EdgeColor', [0.25 0.25 0.25], ...
        'LineWidth', 0.5, ...
        'FaceLighting','none', ...
        'PickableParts','all', ...
        'HitTest','on', ...
        'Tag', sprintf('seg_%d', s));
end
title(ax2, sprintf('Segments: %d   Voxels: %d', nSeg, nnz(voxBin)));

% --- UI: save button & legend hint
uicontrol('Style','pushbutton','String','Save JSON','Units','pixels', ...
          'Position',[10 10 100 28],'Callback',@onSave);
annotation(fig2,'textbox',[0.74 0.01 0.25 0.07], 'String', ...
  sprintf('Hover: highlight\nDouble-click: label\nS: Save   C: Clear\n(Pan/Rotate with mouse)'), ...
  'EdgeColor','none','HorizontalAlignment','right');

% --- Interaction state stored in figure appdata
state.segDir         = segDir;
state.parentDir      = fileparts(segDir);             % one dir ABOVE
state.voxMatPath     = voxMatPath;
state.Hseg           = Hseg;
state.colors         = C;
state.hoverIdx       = 0;
state.selectedMap    = containers.Map('KeyType','int32','ValueType','int32'); % seg->label
state.segSizes       = cellfun(@nnz, segVol);
state.segVoxelCoords = segVoxelCoords;                % for JSON
state.jsonOutPath    = compute_json_path(state.parentDir, voxMatPath); % base_selected_clusters.json
setappdata(fig2, 'state', state);

% --- Callbacks (hover, double-click, keypress)
set(fig2, 'WindowButtonMotionFcn', @onHover);
set(fig2, 'WindowButtonDownFcn',  @onClick);
set(fig2, 'KeyPressFcn',          @onKey);

% ========== nested callbacks ==========
    function onHover(~,~)
        st = getappdata(fig2,'state');
        h  = hittest(fig2);
        k  = hit_to_seg(h, st.Hseg);
        if k == st.hoverIdx, return; end
        % remove old highlight
        if st.hoverIdx>=1 && st.hoverIdx<=numel(st.Hseg) && isvalid(st.Hseg(st.hoverIdx))
            restore_patch_style(st.Hseg(st.hoverIdx), st.colors(st.hoverIdx,:), st.selectedMap);
        end
        % add new highlight
        if k>0
            emphasize_patch(st.Hseg(k));
        end
        st.hoverIdx = k;
        setappdata(fig2,'state',st);
    end

    function onClick(~,~)
        % Label on DOUBLE-CLICK so single-click/drag remains available for pan/rotate
        st = getappdata(fig2,'state');
        if ~strcmp(get(fig2,'SelectionType'),'open') % 'open' = double-click
            return;
        end
        h  = hittest(fig2);
        k  = hit_to_seg(h, st.Hseg);
        if k<=0, return; end
        % prompt small input dialog for integer
        prompt = {sprintf('Label for segment %d (voxels: %d):', k, st.segSizes(k))};
        dlgtitle = 'Assign label';
        dims = [1 35];
        definput = {''};
        answer = inputdlg(prompt, dlgtitle, dims, definput);
        if isempty(answer), return; end
        val = str2double(answer{1});
        if ~isfinite(val) || floor(val)~=val
            warndlg('Please enter an integer.','Invalid input'); return;
        end
        st.selectedMap(int32(k)) = int32(val);
        % show selection by thicker edge
        emphasize_patch(st.Hseg(k), true);
        setappdata(fig2,'state',st);
    end

    function onKey(~,evt)
        st = getappdata(fig2,'state');
        switch lower(evt.Key)
            case 's'   % save
                do_save(st);
            case 'c'   % clear
                st.selectedMap = containers.Map('KeyType','int32','ValueType','int32');
                % restore all styles
                for ii=1:numel(st.Hseg)
                    if isgraphics(st.Hseg(ii))
                        restore_patch_style(st.Hseg(ii), st.colors(ii,:), st.selectedMap);
                    end
                end
                setappdata(fig2,'state',st);
        end
    end

    function onSave(~,~)
        st = getappdata(fig2,'state'); do_save(st);
    end

    function do_save(st)
        % build struct list: {segment_index, voxel_count, label, voxels [x y z]}
        ks = sort(cell2mat(st.selectedMap.keys));
        items = repmat(struct('segment_index',0,'voxel_count',0,'label',0,'voxels',[]), 0, 1);
        for ii = 1:numel(ks)
            k = ks(ii);
            items(end+1) = struct( ... %#ok<AGROW>
                'segment_index', k, ...
                'voxel_count',   st.segSizes(k), ...
                'label',         st.selectedMap(k), ...
                'voxels',        st.segVoxelCoords{k});  % N×3 [x y z]
        end
        out = struct('segments', items, ...
                     'total_segments', numel(st.Hseg), ...
                     'source_segments_dir', st.segDir, ...
                     'source_vox_mat', st.voxMatPath, ...
                     'timestamp', datestr(now,'yyyy-mm-ddTHH:MM:SS'));
        txt = jsonencode(out,'PrettyPrint',true);
    
        % 1) write wherever st.jsonOutPath points (backward-compatible)
        if ~isfield(st,'jsonOutPath') || isempty(st.jsonOutPath)
            st.jsonOutPath = compute_json_path(fileparts(st.segDir), st.voxMatPath); % fallback
        end
        ensure_dir(fileparts(st.jsonOutPath));
        fid = fopen(st.jsonOutPath,'w'); assert(fid>0,'Cannot write JSON: %s', st.jsonOutPath);
        fwrite(fid, txt, 'char'); fclose(fid);
    
        % 2) ensure/copy to canonical annotated_voxel_data path
        canonPath = canonical_json_path(st.voxMatPath);
        if ~strcmp(canonPath, st.jsonOutPath)
            ensure_dir(fileparts(canonPath));
            % prefer move (keeps only one copy). If move fails (e.g., cross-device), fallback to copy.
            [ok,msg] = movefile(st.jsonOutPath, canonPath, 'f');
            if ~ok
                warning('Move failed (%s). Falling back to copy.', msg);
                copyfile(st.jsonOutPath, canonPath, 'f');
            end
            st.jsonOutPath = canonPath; % update state
            setappdata(gcbf,'state',st);
        end
    
        msgbox(sprintf('Saved: %s', canonPath),'Saved','modal');
    end

% ========== helper utils (nested) ==========
    function p = canonical_json_path(voxMatPath)
        % Put JSON in: <viewer_folder>/<base>/annotated_voxel_data/<base>_selected_clusters.json
        here = fileparts(mfilename('fullpath'));          % folder where viewer .m lives
        [~, baseName, ~] = fileparts(voxMatPath);         % e.g., "minW3dmgcg_SMreMF_voxBin"
        suffix = '_voxBin';
        if endsWith(baseName, suffix)
            baseName = extractBefore(baseName, strlength(baseName)-strlength(suffix)+1);
        end
        root_out = fullfile(here, baseName);
        ann_dir  = fullfile(root_out, 'annotated_voxel_data');
        p = fullfile(ann_dir, sprintf('%s_selected_clusters.json', baseName));
    end
    
    function ensure_dir(d)
        if ~isfolder(d), mkdir(d); end
    end

    function k = hit_to_seg(h, H)
        k = 0;
        if ~isscalar(h) || ~isgraphics(h), return; end
        if strcmp(get(h,'Type'),'patch')
            k = find(H==h, 1, 'first');
        end
    end

    function emphasize_patch(ph, isSelected)
        if nargin<2, isSelected=false; end
        if ~isgraphics(ph), return; end
        set(ph, 'EdgeColor', [0 0 0], 'LineWidth', ternary(isSelected,1.8,1.2), 'FaceAlpha', 1.0);
    end

    function restore_patch_style(ph, baseColor, selMap)
        if ~isgraphics(ph), return; end
        idx = hit_to_seg(ph, Hseg);
        isSel = ~isempty(idx) && isKey(selMap, int32(idx));
        set(ph, 'FaceColor', baseColor, ...
                'EdgeColor', [0.25 0.25 0.25], ...
                'LineWidth', ternary(isSel,1.6,0.5), ...
                'FaceAlpha', 1.0);
    end

    function x = ternary(cond,a,b), if cond, x = a; else, x = b; end, end
end % ===== end main =====


% ============================== Local functions ===============================

function voxBin = get_voxbin_from_struct(S)
    if isfield(S,'voxBin')
        voxBin = S.voxBin;
    else
        fns = fieldnames(S);
        ix  = find(contains(lower(fns),'vox'), 1);
        assert(~isempty(ix), 'voxBin variable not found in MAT file.');
        voxBin = S.(fns{ix});
    end
    assert(ndims(voxBin)==3, 'voxBin must be 3-D.');
    voxBin = logical(voxBin);
end

function V = load_obj_vertices(path)
    fid = fopen(path,'r');
    assert(fid>0, 'Failed to open OBJ: %s', path);
    C = onCleanup(@() fclose(fid));
    verts = [];
    while true
        L = fgetl(fid);
        if ~ischar(L), break; end
        if numel(L)>=2 && (L(1)=='v' && isspace(L(2)))
            nums = sscanf(L(2:end), '%f');
            if numel(nums)>=3
                verts(end+1,1:3) = nums(1:3).'; %#ok<AGROW>
            end
        end
    end
    if isempty(verts), V = zeros(0,3); else, V = verts; end
end

function [segVol, vox2seg] = map_voxels_to_segments(voxBin, voxCoord, Vseg)
% True line-segment distance from each voxel center to each polyline edge (chunked)
    nSeg = numel(Vseg);
    % Build edge list
    segA = []; segB = []; segIDofEdge = [];
    for s = 1:nSeg
        V = Vseg{s};
        if size(V,1) >= 2
            segA = [segA; V(1:end-1,:)]; %#ok<AGROW>
            segB = [segB; V(2:end  ,:)]; %#ok<AGROW>
            segIDofEdge = [segIDofEdge; repmat(s, size(V,1)-1, 1)]; %#ok<AGROW>
        end
    end
    E   = size(segA,1);
    BA  = segB - segA;                   % (E×3)
    BA2 = sum(BA.^2, 2).';               % (1×E)

    Nvox = size(voxCoord,1);
    vox2seg = zeros(Nvox, 1, 'uint32');

    block = 5e4;
    for i = 1:block:Nvox
        j  = min(i+block-1, Nvox);
        P  = double(voxCoord(i:j,:));    % (b×3)
        b  = size(P,1);
        P3  = reshape(P,  [b 1 3]);
        A3  = reshape(segA,[1 E 3]);
        BA3 = reshape(BA, [1 E 3]);
        t = sum( (P3 - A3) .* BA3, 3) ./ BA2;
        t = max(0, min(1, t));
        proj = A3 + BA3 .* reshape(t,[b E 1]);
        d2   = sum( (P3 - proj).^2, 3 );
        [~,edgeIdx] = min(d2, [], 2);
        vox2seg(i:j) = uint32(segIDofEdge(edgeIdx));
    end

    % Convert mapping to per-segment masks
    [idxY, idxZ, idxX] = ind2sub(size(voxBin), find(voxBin));
    segVol = cell(1, nSeg);
    for s = 1:nSeg
        mask = false(size(voxBin));
        sel  = (vox2seg == s);
        if any(sel)
            lin = sub2ind(size(voxBin), idxY(sel), idxZ(sel), idxX(sel));
            mask(lin) = true;
        end
        segVol{s} = mask;
    end
end

function [R, t] = choose_alignment(voxBin, Vseg)
% AutoAlign v3: weighted centroid + overlap scoring + tiny integer refine
    [ny,nz,nx] = size(voxBin);
    voxMask = voxBin;
    try, voxMask = imdilate(voxBin, true(3)); catch, end
    [yy,zz,xx] = ind2sub(size(voxBin), find(voxBin));
    voxXYZ = [xx(:) yy(:) zz(:)];
    c_vox = mean(voxXYZ,1);
    try
        D = bwdist(~voxMask); w = double(D(voxMask));
        if ~isempty(w) && any(w>0)
            c_vox = sum(voxXYZ.*w,1) ./ sum(w);
        end
    catch, end
    sk_all = vertcat(Vseg{:});
    Rz = cat(3, eye(3), [0 -1 0; 1 0 0; 0 0 1], [-1 0 0; 0 -1 0; 0 0 1], [0 1 0; -1 0 0; 0 0 1]);
    bestScore = -inf; bestR = eye(3); bestt = [0 0 0];
    function sc = overlap_score(W2)
        idx = round(W2);
        idx(:,1) = min(max(idx(:,1),1), nx);
        idx(:,2) = min(max(idx(:,2),1), ny);
        idx(:,3) = min(max(idx(:,3),1), nz);
        lin = sub2ind([ny,nz,nx], idx(:,2), idx(:,3), idx(:,1));
        sc = mean(voxMask(lin));
    end
    for i = 1:size(Rz,3)
        R0 = Rz(:,:,i);
        W  = (R0 * sk_all.').';
        c_W = mean(W,1);
        t0  = c_vox - c_W;
        bestLocalScore = -inf; bestLocalT = t0;
        for dx = -5:5
            for dy = -1:1
                for dz = -1:1
                    tv = t0 + [dx dy dz];
                    sc = overlap_score(bsxfun(@plus, W, tv));
                    if sc > bestLocalScore
                        bestLocalScore = sc; bestLocalT = tv;
                    end
                end
            end
        end
        if bestLocalScore > bestScore
            bestScore = bestLocalScore; bestR = R0; bestt = bestLocalT;
        end
    end
    R = bestR; t = bestt;
end

function C = simple_distinguishable_colors(n)
    h = linspace(0,1,n+1); h(end) = []; s = 0.75; v = 0.95;
    C = hsv2rgb([h(:), repmat(s,n,1), repmat(v,n,1)]);
    C = C(mod(round((0:n-1)*2.3), n)+1, :);
end

function [V,F] = build_voxel_mesh_quads(mask)
% Quad mesh of exposed voxel faces for binary volume "mask" (ny×nz×nx, y-z-x).
    [ny,nz,nx] = size(mask);
    xp = false(ny,nz,nx); xp(:,:,1:end-1) = ~mask(:,:,2:end) & mask(:,:,1:end-1);
    xm = false(ny,nz,nx); xm(:,:,2:end)   = ~mask(:,:,1:end-1) & mask(:,:,2:end);
    yp = false(ny,nz,nx); yp(1:end-1,:,:) = ~mask(2:end,:,:)   & mask(1:end-1,:,:);
    ym = false(ny,nz,nx); ym(2:end,:,:)   = ~mask(1:end-1,:,:) & mask(2:end,:,:);
    zp = false(ny,nz,nx); zp(:,1:end-1,:) = ~mask(:,2:end,:)   & mask(:,1:end-1,:);
    zm = false(ny,nz,nx); zm(:,2:end,:)   = ~mask(:,1:end-1,:) & mask(:,2:end,:);

    V = []; F = [];
    function add_face(quadVerts)
        base = size(V,1);
        V = [V; quadVerts]; %#ok<AGROW>
        F = [F; base+(1:4)]; %#ok<AGROW>
    end

    [yy,zz,xx] = ind2sub(size(xp), find(xp));
    for k = 1:numel(xx)  % +X face at x+0.5
        x = xx(k); y = yy(k); z = zz(k);
        x0 = x+0.5; y0 = y; z0 = z;
        add_face([ x0 y0-0.5 z0-0.5;  x0 y0+0.5 z0-0.5;  x0 y0+0.5 z0+0.5;  x0 y0-0.5 z0+0.5 ]);
    end
    [yy,zz,xx] = ind2sub(size(xm), find(xm));
    for k = 1:numel(xx)  % -X face
        x = xx(k); y = yy(k); z = zz(k);
        x0 = x-0.5; y0 = y; z0 = z;
        add_face([ x0 y0-0.5 z0+0.5;  x0 y0+0.5 z0+0.5;  x0 y0+0.5 z0-0.5;  x0 y0-0.5 z0-0.5 ]);
    end
    [yy,zz,xx] = ind2sub(size(yp), find(yp));
    for k = 1:numel(xx)  % +Y face
        x = xx(k); y = yy(k); z = zz(k);
        y0 = y+0.5; x0 = x; z0 = z;
        add_face([ x0-0.5 y0 z0-0.5;  x0+0.5 y0 z0-0.5;  x0+0.5 y0 z0+0.5;  x0-0.5 y0 z0+0.5 ]);
    end
    [yy,zz,xx] = ind2sub(size(ym), find(ym));
    for k = 1:numel(xx)  % -Y face
        x = xx(k); y = yy(k); z = zz(k);
        y0 = y-0.5; x0 = x; z0 = z;
        add_face([ x0-0.5 y0 z0+0.5;  x0+0.5 y0 z0+0.5;  x0+0.5 y0 z0-0.5;  x0-0.5 y0 z0-0.5 ]);
    end
    [yy,zz,xx] = ind2sub(size(zp), find(zp));
    for k = 1:numel(xx)  % +Z face
        x = xx(k); y = yy(k); z = zz(k);
        z0 = z+0.5; x0 = x; y0 = y;
        add_face([ x0-0.5 y0-0.5 z0;  x0+0.5 y0-0.5 z0;  x0+0.5 y0+0.5 z0;  x0-0.5 y0+0.5 z0 ]);
    end
    [yy,zz,xx] = ind2sub(size(zm), find(zm));
    for k = 1:numel(xx)  % -Z face
        x = xx(k); y = yy(k); z = zz(k);
        z0 = z-0.5; x0 = x; y0 = y;
        add_face([ x0-0.5 y0+0.5 z0;  x0+0.5 y0+0.5 z0;  x0+0.5 y0-0.5 z0;  x0-0.5 y0-0.5 z0 ]);
    end
end

function outPath = compute_json_path(parentDir, voxMatPath)
% Save JSON one directory ABOVE segments folder, name it "<base>_selected_clusters.json"
    [~, baseName, ~] = fileparts(voxMatPath);  % e.g., "minW3dmgcg_SMreMF_voxBin"
    suffix = '_voxBin';
    if endsWith(baseName, suffix), baseName = extractBefore(baseName, strlength(baseName)-strlength(suffix)+1); end
    outPath = fullfile(parentDir, sprintf('%s_selected_clusters.json', baseName));
end
