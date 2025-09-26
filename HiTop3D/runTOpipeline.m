function runTOpipeline(nelx,nely,nelz,penal,rmin,ft,nl,cgtol,cgmax,compconst,filename, varargin)
% runTOpipeline — TO → STL → cap → separate (LS skeleton) → segment → verify → organize
%
% Usage:
%   >> runTOpipeline(48,24,24,3,sqrt(3),1,4,1e0,100,3330,'minW3dmgcg_SMreMF')
%
% Options:
%   'SkipExisting' (logical, default true)
%   'Force'        (logical, default false)
%   'Visualize'    (logical, default true)
%   'Python'       (char,    default auto-detect)
%   'AutoPip'      (logical, default false)
%   'PipArgs'      (char,    default '')

    % ------------------ options ------------------
    opts = struct('SkipExisting', true, 'Force', false, 'Visualize', true, ...
                  'Python','', 'AutoPip', false, 'PipArgs','');
    opts = parseOptions(opts, varargin{:});
    if opts.Force, opts.SkipExisting = false; end

    PYTHON_CMD = opts.Python;
    if isempty(PYTHON_CMD), PYTHON_CMD = detect_python(); end

    here = fileparts(mfilename('fullpath'));
    addpath(here);

    % ------------------ script paths ------------------
    path_minW      = fullfile(here,'minW3dmgcg_SMreMF.m');
    path_vox2stl   = fullfile(here,'voxBin2stl.m');
    path_view      = fullfile(here,'viewSkeletonVsVoxel.m');
    py_quickcap    = fullfile(here,'quick_cap.py');
    py_separator   = fullfile(here,'separator.py');  % Step 4

    % ------------------ names / I/O (original working-dir paths) ----------
    base      = filename;
    voxMat    = [base, '_voxBin.mat'];                  % Step 1
    stl_raw   = [base, '.stl'];                         % Step 2
    stl_cap   = [base, '_capped.stl'];                  % Step 3
    skel_obj  = [base, '_capped_skel_ls.obj'];          % Step 4 output
    seg_dir   = [base, '_capped_skel_ls_segments'];     % Step 5

    % ------------------ NEW: organized roots & organized paths ------------
    root_out = fullfile(here, base);                    % keep outputs next to scripts
    ann_dir  = fullfile(root_out, 'annotated_voxel_data');
    skel_out = fullfile(root_out, 'skeletonization_data');
    ensure_dir(root_out); ensure_dir(ann_dir); ensure_dir(skel_out);

    voxMat_ann    = fullfile(ann_dir, [base, '_voxBin.mat']);
    stl_raw_skel  = fullfile(skel_out, [base, '.stl']);
    stl_cap_skel  = fullfile(skel_out, [base, '_capped.stl']);
    skel_obj_skel = fullfile(skel_out, [base, '_capped_skel_ls.obj']);
    seg_dir_skel  = fullfile(skel_out, [base, '_capped_skel_ls_segments']);

    % ------------------ NEW: resolve where artifacts live now -------------
    voxMat   = resolve_artifact(voxMat,   voxMat_ann);
    stl_raw  = resolve_artifact(stl_raw,  stl_raw_skel);
    stl_cap  = resolve_artifact(stl_cap,  stl_cap_skel);
    skel_obj = resolve_artifact(skel_obj, skel_obj_skel);
    if isfolder(seg_dir_skel), seg_dir = seg_dir_skel; end

    % JSON path the viewer may create (organizer will move it if not already here)
    selected_json = fullfile(ann_dir, sprintf('%s_selected_clusters.json', base));

    % For step headings
    N_STEPS = 7;

    % ------------------ Step 0: Python deps preflight --
    stepHeader(0,N_STEPS, 'Python dependency check for quick_cap.py');
    python_info(PYTHON_CMD);
    ensure_py_modules(PYTHON_CMD, {'meshio','pymeshfix'}, opts);

    % ------------------ Step 1: TO --------------------
    stepHeader(1,N_STEPS, sprintf('minW3dmgcg_SMreMF → %s', voxMat));
    must_exist(path_minW, 'minW3dmgcg_SMreMF.m not found.');
    if shouldSkipAny(opts, {voxMat, voxMat_ann})      % NEW
        fprintf('  • Skip (exists)\n');
    else
        minW3dmgcg_SMreMF(nelx,nely,nelz,penal,rmin,ft,nl,cgtol,cgmax,compconst,base);
        must_exist([base, '_voxBin.mat'], 'voxBin output missing after TO step.');
        voxMat = [base, '_voxBin.mat']; % created in working dir; organizer may move later
    end

    % ------------------ Step 2: vox → STL -------------
    stepHeader(2,N_STEPS, sprintf('voxBin2stl → %s', stl_raw));
    must_exist(path_vox2stl, 'voxBin2stl.m not found.');
    must_exist(voxMat,       'voxBin input missing.');
    if shouldSkipAny(opts, {stl_raw, stl_raw_skel})    % NEW
        fprintf('  • Skip (exists)\n');
    else
        voxBin2stl(voxMat, [base, '.stl']);
        must_exist([base, '.stl'], 'STL missing after voxBin2stl.');
        stl_raw = [base, '.stl'];
    end

    % ------------------ Step 3: cap STL ----------------
    stepHeader(3,N_STEPS, sprintf('quick_cap.py → %s', stl_cap));
    must_exist(py_quickcap, 'quick_cap.py missing.');
    must_exist(stl_raw,     'raw STL input missing.');
    if shouldSkipAny(opts, {stl_cap, stl_cap_skel})    % NEW
        fprintf('  • Skip (exists)\n');
    else
        ensure_py_modules(PYTHON_CMD, {'meshio','pymeshfix'}, opts);
        cmd_cap = sprintf('"%s" "%s" "%s" "%s"', PYTHON_CMD, py_quickcap, stl_raw, [base, '_capped.stl']);
        run_sys(cmd_cap);
        must_exist([base, '_capped.stl'], 'Capped STL missing.');
        stl_cap = [base, '_capped.stl'];
    end

    % ------------------ Step 4: separator (LS skeleton) ----
    stepHeader(4,N_STEPS, sprintf('separator.py → %s', skel_obj));
    must_exist(py_separator, 'separator.py missing.');
    must_exist(stl_cap,      'capped STL input missing.');

    ensure_py_modules(PYTHON_CMD, {'numpy','trimesh','pygel3d'}, opts);

    if shouldSkipAny(opts, {skel_obj, skel_obj_skel})  % NEW
        fprintf('  • Skip (exists)\n');
    else
        wrapper_sep = fullfile(tempdir, ['sep_wrap_' char(java.util.UUID.randomUUID) '.py']);
        write_text(wrapper_sep, build_separator_wrapper(py_separator, stl_cap, [base, '_capped_skel_ls.obj']));
        cmd_sep = sprintf('"%s" "%s"', PYTHON_CMD, wrapper_sep);
        run_sys(cmd_sep);
        must_exist([base, '_capped_skel_ls.obj'], 'Skeleton OBJ missing after separator.py.');
        skel_obj = [base, '_capped_skel_ls.obj'];
    end

    % ------------------ Step 5: segmentor -------------
    stepHeader(5,N_STEPS, sprintf('segmentor.py → %s/', seg_dir));
    py_segmentor = fullfile(here,'segmentor.py');
    must_exist(py_segmentor, 'segmentor.py missing.');
    must_exist(skel_obj,     'skeleton OBJ input missing.');
    if shouldSkipAny(opts, {seg_dir, seg_dir_skel})    % NEW
        fprintf('  • Skip (exists)\n');
    else
        wrapper2 = fullfile(tempdir, ['seg_wrap_' char(java.util.UUID.randomUUID) '.py']);
        write_text(wrapper2, build_segmentor_wrapper(py_segmentor, skel_obj));
        cmd_seg = sprintf('"%s" "%s"', PYTHON_CMD, wrapper2);
        run_sys(cmd_seg);
        if ~isfolder([base, '_capped_skel_ls_segments'])
            error('Segments folder not found: %s', [base, '_capped_skel_ls_segments']);
        end
        seg_dir = [base, '_capped_skel_ls_segments'];
    end

    % ------------------ Step 6: visualize -------------
    stepHeader(6,N_STEPS, sprintf('viewSkeletonVsVoxel("%s", "%s")', seg_dir, voxMat));
    if opts.Visualize
        must_exist(path_view, 'viewSkeletonVsVoxel.m not found.');
        must_exist(voxMat,    'voxBin input missing.');
        if ~isfolder(seg_dir)
            warning('Segments folder not found; skipping visualization: %s', seg_dir);
        else
            % Viewer will save JSON; our organizer will ensure it ends up under ann_dir.
            viewSkeletonVsVoxel(seg_dir, voxMat);
        end
    else
        fprintf('  • Skip (disabled)\n');
    end

    % ------------------ Step 7: organize outputs ------
    stepHeader(7,N_STEPS, 'Organize outputs (annotated_voxel_data/ & skeletonization_data/)');

    ensure_dir(root_out);
    ensure_dir(ann_dir);
    ensure_dir(skel_out);

    % 7a) Annotated voxel data: move all MATs for this base + selected_clusters.json
    mats = dir([base, '*.mat']);
    for k = 1:numel(mats)
        move_if_exists(fullfile(mats(k).folder, mats(k).name), fullfile(ann_dir, mats(k).name), opts);
    end
    % If viewer saved JSON elsewhere, ensure it's in annotated_voxel_data
    other_json = fullfile(fileparts(seg_dir), sprintf('%s_selected_clusters.json', base));
    if isfile(other_json) && ~strcmp(other_json, selected_json)
        move_if_exists(other_json, selected_json, opts);
    end
    if isfile(selected_json)
        % already correct place (no-op), but ensures path exists for downstream
    end

    % 7b) Skeletonization data: STL/OBJ/segments + any other skeleton JSONs (except selected_clusters.json)
    move_if_exists(stl_raw,  fullfile(skel_out, getname(stl_raw)),  opts);
    move_if_exists(stl_cap,  fullfile(skel_out, getname(stl_cap)),  opts);
    move_if_exists(skel_obj, fullfile(skel_out, getname(skel_obj)), opts);

    if isfolder(seg_dir)
        move_dir_if_exists(seg_dir, fullfile(skel_out, getname(seg_dir)), opts);
        seg_dir = fullfile(skel_out, getname(seg_dir)); % update resolved location
    end

    % Move other JSONs (exclude selected_clusters.json already handled)
    js = dir('*.json');
    for k = 1:numel(js)
        src = fullfile(js(k).folder, js(k).name);
        if strcmpi(src, selected_json), continue; end
        move_if_exists(src, fullfile(skel_out, js(k).name), opts);
    end

    % ------------------ NEW: final re-resolve for future runs -------------
    voxMat   = resolve_artifact(voxMat,   voxMat_ann);
    stl_raw  = resolve_artifact(stl_raw,  stl_raw_skel);
    stl_cap  = resolve_artifact(stl_cap,  stl_cap_skel);
    skel_obj = resolve_artifact(skel_obj, skel_obj_skel);
    if isfolder(seg_dir_skel), seg_dir = seg_dir_skel; end

    fprintf('\n✓ Pipeline completed.\n');
    fprintf('  • Outputs organized under: %s\n', root_out);
    fprintf('    - %s\n', ann_dir);
    fprintf('    - %s\n', skel_out);
end

% ---------- helpers ----------
function yes = shouldSkipAny(opts, paths)
    % paths: cellstr of files or folders; true if any exists and skipping enabled
    if ~opts.SkipExisting
        yes = false; return;
    end
    yes = false;
    for i = 1:numel(paths)
        p = paths{i};
        if isfile(p) || isfolder(p)
            yes = true; return;
        end
    end
end

function p = resolve_artifact(primaryPath, organizedPath)
    % Prefer whichever exists; otherwise default to primaryPath (where the step will write)
    if isfile(organizedPath) || isfolder(organizedPath)
        p = organizedPath;
    elseif isfile(primaryPath) || isfolder(primaryPath)
        p = primaryPath;
    else
        p = primaryPath; % not created yet
    end
end

function stepHeader(i,n,msg)
    if i==0
        fprintf('\n[Preflight] %s\n', msg);
    else
        fprintf('\n[%d/%d] %s\n', i, n, msg);
    end
end

function must_exist(p,msg)
    if ~isfile(p) && ~isfolder(p), error('%s  [%s]', msg, p); end
end

function run_sys(cmd)
    [st, out] = system(cmd);
    if st ~= 0
        fprintf(2, '\nCommand failed:\n%s\n', cmd);
        if ~isempty(out), fprintf(2, '%s\n', out); end
        error('External command returned error.');
    else
        if ~isempty(out), fprintf('%s\n', out); end
    end
end

function write_text(path, txt)
    fid = fopen(path, 'w');  assert(fid>0, 'Cannot create wrapper: %s', path);
    cleaner = onCleanup(@() fclose(fid));
    fwrite(fid, txt, 'char');
end

function exe = detect_python()
    if ispc
        exe = 'python';
    else
        [st,~] = system('command -v python3');
        if st==0, exe = 'python3'; else, exe = 'python'; end
    end
end

function r = ternary(cond,a,b)
    if cond, r=a; else, r=b; end
end

function opts = parseOptions(opts, varargin)
    if mod(numel(varargin),2)~=0
        error('Options must be name-value pairs.');
    end
    for k=1:2:numel(varargin)
        name = varargin{k}; val = varargin{k+1};
        ln = lower(char(string(name)));
        switch ln
            case 'skipexisting', opts.SkipExisting = logical(val);
            case 'force',        opts.Force       = logical(val);
            case 'visualize',    opts.Visualize   = logical(val);
            case 'python',       opts.Python      = char(val);
            case 'autopip',      opts.AutoPip     = logical(val);
            case 'pipargs',      opts.PipArgs     = char(val);
            otherwise, error('Unknown option: %s', name);
        end
    end
end

% ===== Robust JSON quoting for Python wrappers =====
function j = json_quote(s)
    s = strrep(s, '\', '\\');
    s = strrep(s, '"', '\"');
    j = ['"', s, '"'];
end

function code = build_separator_wrapper(py_separator, stl_in, expected_skel_out)
    code = sprintf([ ...
        'import importlib.util, pathlib\n' ...
        'p = %s\n' ...
        'spec = importlib.util.spec_from_file_location("sep", p)\n' ...
        'm = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)\n' ...
        'm.PATH_MESH = %s\n' ...
        'm.BASENAME  = pathlib.Path(m.PATH_MESH).with_suffix("""""").name + "_skel"\n' ...
        'm.main()\n' ...
        'out_expected = %s\n' ...
        'from pathlib import Path as _P\n' ...
        'assert _P(out_expected).exists(), f"Expected skeleton not found: {out_expected}"\n' ...
        ], json_quote(py_separator), json_quote(stl_in), json_quote(expected_skel_out));
end

function code = build_segmentor_wrapper(py_segmentor, skel_obj)
    code = sprintf([ ...
        'import importlib.util\n' ...
        'p = %s\n' ...
        'spec = importlib.util.spec_from_file_location("segm", p)\n' ...
        'm = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)\n' ...
        'm.INPUT_OBJ = %s\n' ...
        'm.main()\n' ...
        ], json_quote(py_segmentor), json_quote(skel_obj));
end

% ===== Python preflight & info =====
function ensure_py_modules(pythonExe, mods, opts)
    % mods: cellstr of pip package names
    safeCwd = tempdir;

    function name = import_name_for(pipname)
        switch lower(pipname)
            case 'scikit-image',  name = 'skimage';
            case 'opencv-python', name = 'cv2';
            case 'pyyaml',        name = 'yaml';
            otherwise,            name = regexprep(pipname,'-','_');
        end
    end

    function st = can_import(modname)
        checkCmd = sprintf( ...
            'cd "%s"; "%s" -c "ok=1\ntry:\n import %s\nexcept Exception:\n ok=0\nimport sys; sys.exit(0 if ok else 1)"', ...
            safeCwd, pythonExe, modname);
        st = system(checkCmd);
    end

    missing_pip = {};
    for i = 1:numel(mods)
        pipname = mods{i};
        impname = import_name_for(pipname);
        if can_import(impname) ~= 0
            missing_pip{end+1} = pipname; %#ok<AGROW>
        end
    end

    if isempty(missing_pip)
        fprintf('  • Python deps OK (%s)\n', strjoin(mods, ', '));
        return;
    end

    fprintf(2, '  • Missing Python modules: %s\n', strjoin(missing_pip, ', '));
    if ~opts.AutoPip
        fprintf(2, '    Re-run with ''AutoPip'', true or install manually for THIS Python:\n');
        fprintf(2, '      %s -m pip install %s %s\n', pythonExe, strjoin(missing_pip, ' '), opts.PipArgs);
        error('Missing required Python modules.');
    end

    pipCmd = sprintf('"%s" -m pip install %s %s', pythonExe, strjoin(missing_pip, ' '), opts.PipArgs);
    fprintf('  • Installing with: %s\n', pipCmd);
    [st, out] = system(pipCmd);
    if st ~= 0
        fprintf(2, '%s\n', out);
        error('AutoPip failed. Try installing manually.');
    end

    for i = 1:numel(missing_pip)
        impname = import_name_for(missing_pip{i});
        if can_import(impname) ~= 0
            error('Module "%s" still missing after pip install (import "%s" failed).', missing_pip{i}, impname);
        end
    end
    fprintf('  • Install OK.\n');
end

function python_info(pythonExe)
    fprintf('  • Python exe / version / site:\n');
    info_py = { ...
        'import sys, site', ...
        'print(sys.executable)', ...
        'print(sys.version.replace(chr(10), '' ''))', ...
        'print(''site-packages:'')', ...
        'paths = []', ...
        'try:', ...
        '    paths += site.getsitepackages()', ...
        'except Exception:', ...
        '    pass', ...
        'try:', ...
        '    paths.append(site.getusersitepackages())', ...
        'except Exception:', ...
        '    pass', ...
        'print(''\n''.join(paths))' ...
    };
    script = strjoin(info_py, sprintf('\n'));
    tmp = fullfile(tempdir, ['pyinfo_' char(java.util.UUID.randomUUID) '.py']);
    write_text(tmp, script);
    cmd = sprintf('"%s" "%s"', pythonExe, tmp);
    [~, out] = system(cmd);
    fprintf('%s\n', out);
end

% ===== filesystem helpers =====
function ensure_dir(d)
    if ~isfolder(d), mkdir(d); end
end

function move_if_exists(src, dst, opts)
    if ~isfile(src), return; end
    dstDir = fileparts(dst);
    if ~isfolder(dstDir), mkdir(dstDir); end
    if isfile(dst)
        if opts.Force
            delete(dst);
        else
            % keep existing — write alongside with suffix
            [p,n,e] = fileparts(dst);
            dst = fullfile(p, sprintf('%s_existing%s', n, e));
        end
    end
    [ok,msg] = movefile(src, dst, 'f');
    if ~ok, warning('Failed to move "%s" → "%s": %s', src, dst, msg); end
end

function move_dir_if_exists(srcDir, dstDir, opts)
    if ~isfolder(srcDir), return; end
    if isfolder(dstDir)
        if opts.Force
            rmdir(dstDir, 's');
        else
            % write alongside with suffix
            dstDir = [dstDir, '_existing'];
        end
    end
    parent = fileparts(dstDir);
    if ~isfolder(parent), mkdir(parent); end
    [ok,msg] = movefile(srcDir, dstDir);
    if ~ok, warning('Failed to move folder "%s" → "%s": %s', srcDir, dstDir, msg); end
end

function n = getname(p)
    [~,n,ext] = fileparts(p);
    n = [n, ext];
end

