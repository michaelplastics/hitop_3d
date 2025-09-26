%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% voxBin2stl.m  —  Voxels → watertight surface mesh  (binary STL)
%
%   voxBin2stl('model_voxBin.mat');              % → model.stl
%   voxBin2stl('model_voxBin.mat','mesh.stl');   % custom name
%
% REQUIREMENTS
%   • Image Processing Toolbox     (bwdist, imfill, smooth3, isosurface)
%   • stlwrite (built-in ≥ R2018a)  – or File-Exchange fallback
%   • Iso2Mesh optional: meshcheckrepair  (for extra watertightness)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function voxBin2stl(voxMat, stlOut)

% ---------- filenames ----------------------------------------------------
if nargin<1
    error('Usage: voxBin2stl(''file_voxBin.mat'' [, ''mesh.stl''])');
end
if nargin<2
    [p,b,~] = fileparts(voxMat);
    b = regexprep(b,'_voxBin$','');
    stlOut = fullfile(p,[b '.stl']);
end

% ---------- load voxels --------------------------------------------------
S = load(voxMat);
assert(isfield(S,'voxBin'),'Variable ''voxBin'' not found in %s',voxMat);
voxBin = logical(S.voxBin);         % [y z x]

% ---------- seal holes in volume ----------------------------------------
vol = permute(voxBin,[3 1 2]);      % → x y z
vol = padarray(vol,[1 1 1],0);      % pad so fill can escape
vol = imfill(vol,'holes');          % flood-fill 3-D (IPT R2018a+)
vol = vol(2:end-1,2:end-1,2:end-1); % remove padding

% ---------- optional slight smoothing -----------------------------------
%vol = smooth3(vol,'box',3);

% ---------- extract surface ---------------------------------------------
isoVal = 0.5;
[F,V] = isosurface(vol, isoVal);     % triangulated faces / verts

% ---------- optional watertight repair ----------------------------------
if exist('meshcheckrepair','file') == 2       % Iso2Mesh present
    [V,F] = meshcheckrepair(V,F,'dup');       % remove dups & fill holes
end

% ---------- write binary STL --------------------------------------------
fprintf('Writing %s  (verts:%d  faces:%d)\n', stlOut, size(V,1), size(F,1));
try
    % newer syntax (R2019b +)
    TR = triangulation(F, V);
    stlwrite(TR, stlOut);          % stlwrite(TR,'file.stl')
catch
    % older syntax
    stlwrite(stlOut, F, V);        % stlwrite('file.stl',F,V)
end

fprintf('Done.\n');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
