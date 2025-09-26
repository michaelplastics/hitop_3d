function mapVoxelsToSkeleton(objFolder, voxMatFile, outPrefix)
% mapVoxelsToSkeleton
% -------------------------------------------------------------------------
%  • Reads every *.obj in objFolder   (one file = one skeleton segment)
%  • Loads voxBin (y,z,x logical) from voxMatFile
%  • Assigns each solid voxel to the nearest segment (nearest OBJ vertex)
%  • Saves:
%        <outPrefix>_segments_voxBin.mat   % cell array of logical volumes
%        <outPrefix>_segXX_vox.mat         % one file per segment
%  • Opens a 3-D scatter preview – each segment in a different colour
%
% USAGE
%   mapVoxelsToSkeleton('Segments', 'voxBin.mat');           % default prefix
%   mapVoxelsToSkeleton('Segments', 'voxBin.mat', 'mapped'); % custom prefix
%
% REQUIREMENTS
%   – Pure MATLAB (no toolboxes).  Uses pdist2 / knnsearch if Statistics
%     Toolbox is present, otherwise a small loop-based fallback.
% -------------------------------------------------------------------------
if nargin < 3
    [~,b,~] = fileparts(voxMatFile);
    outPrefix = b;
end
assert(exist(objFolder,'dir')==7,  'OBJ folder not found.');
assert(exist(voxMatFile,'file')==2,'voxMatFile not found.');

%% 1.  read skeleton segments (OBJ vertices) ------------------------------
objFiles = dir(fullfile(objFolder,'*.obj'));
assert(~isempty(objFiles),'No *.obj files in %s',objFolder);

% rigid transform that aligns the skeleton to voxBin
R = [ 0  1  0;           % 90° clockwise about +Z
     -1  0  0;
      0  0  1];
t = [0  24  0];          % +24 voxels along +Y

Vseg   = cell(1,numel(objFiles));     % each cell = [n_i × 3] vertices
segLab = [];                          % segment-label per vertex
Vall   = [];                          % stacked vertices

for s = 1:numel(objFiles)
    V = readOBJvertices(fullfile(objFiles(s).folder,objFiles(s).name));
    V = (R * V.').' + t;              % apply rotation AND translation
    Vseg{s} = V;                      %#ok<*AGROW>
    Vall    = [Vall; V];
    segLab  = [segLab; repmat(s,size(V,1),1)];
end
nSeg = numel(Vseg);
fprintf('Loaded %d segments  (total vertices: %d)\n', ...
        nSeg, size(Vall,1));

%% 2.  load voxels ---------------------------------------------------------
S = load(voxMatFile);  assert(isfield(S,'voxBin'),'voxBin not in MAT');
voxBin = logical(S.voxBin);                 % [y z x]
[ny,nz,nx] = size(voxBin);

[idxY,idxZ,idxX] = ind2sub(size(voxBin), find(voxBin));
voxCoord = [idxX idxY idxZ];                % (N × 3) in XYZ

% ---------- toolbox-free closest-segment mapping -------------------------
block = 5e4;                             % voxels per chunk (adjust to RAM)
vox2seg = zeros(size(voxCoord,1),1,'uint32');

segA = Vseg{1}(1,:);  % dummy init
segA = []; segB = [];
for s = 1:nSeg
    V   = Vseg{s};
    segA = [segA; V(1:end-1,:)];         % start of every polyline edge
    segB = [segB; V(2:end  ,:)];         % end   of every edge
end
BA   = segB - segA;                      % (E×3)
BA2  = sum(BA.^2,2).';                   % 1×E   (squared lengths)
segIDofEdge = repelem(1:nSeg, ...
                      cellfun(@(v)size(v,1)-1,Vseg)).';

E = size(segA,1);                        % # edges

for i = 1:block:size(voxCoord,1)
    j  = min(i+block-1, size(voxCoord,1));
    P  = double(voxCoord(i:j,:));        % (b×3)
    b  = size(P,1);

    % reshape to (b × E × 3)
    P3   = reshape(P,  [b 1 3]);
    A3   = reshape(segA, [1 E 3]);
    BA3  = reshape(BA,  [1 E 3]);

    % projection parameter t in 0…1 along each segment
    t = sum( (P3 - A3) .* BA3, 3) ./ BA2;     % (b×E)
    t = max(0, min(1, t));

    proj = A3 + BA3 .* reshape(t,[b E 1]);     % (b×E×3)
    d2   = sum( (P3 - proj).^2, 3 );           % (b×E)

    [~,edgeIdx] = min(d2,[],2);                % closest edge per voxel
    vox2seg(i:j) = segIDofEdge(edgeIdx);
end
fprintf('Voxel→segment mapping done (true line distance).\n');

% ---------- 4. build logical volumes per segment ------------------------
segVol = cell(1,nSeg);

segDir = [outPrefix '_segments'];        % <-- NEW
if ~exist(segDir,'dir'), mkdir(segDir); end

for s = 1:nSeg
    mask = false(size(voxBin));
    mask( sub2ind([ny nz nx], ...
        idxY(vox2seg==s), idxZ(vox2seg==s), idxX(vox2seg==s)) ) = true;
    segVol{s} = mask;

    save(fullfile(segDir, sprintf('seg%02d_vox.mat',s)), ...  % <-- NEW
         'mask');
end
save(fullfile(segDir,[outPrefix '_segments_voxBin.mat']), ... % <-- NEW
     'segVol','-v7');
fprintf('Saved %d segment volumes to folder %s\n', nSeg, segDir);
% --- sanity report -------------------------------------------------------
voxPerSeg = cellfun(@nnz, segVol);
fprintf('\nVoxel count per segment:\n');
disp([ (1:nSeg).' voxPerSeg.' ])          % two-column list

emptySeg  = find(voxPerSeg==0);
if ~isempty(emptySeg)
    warning('Segments with ZERO voxels: %s\n', num2str(emptySeg));
end
%% 5.  preview with large "voxel" squares ---------------------------------
C = distinguishable_colors(nSeg, [1 1 1]);
mk = 150;                     % <-- marker size in pixels (try 100–300)
figure('Name','Segment mapping'); hold on
for s = 1:nSeg
    [y,z,x] = ind2sub(size(voxBin), find(segVol{s}));
    scatter3(x, y, z, mk, C(s,:), 's', 'filled');   % 's' = square marker
end
axis equal ij off; view(3);
title(sprintf('Segments: %d   Total voxels: %d', nSeg, nnz(voxBin)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V = readOBJvertices(filename)
% very light OBJ vertex parser (ignores faces)
fid = fopen(filename,'r');
assert(fid>0,'Cannot open %s',filename);
V = textscan(fid,'v %f %f %f%*[^\n]','CollectOutput',true, ...
             'CommentStyle','#');
fclose(fid);
V = V{1};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%