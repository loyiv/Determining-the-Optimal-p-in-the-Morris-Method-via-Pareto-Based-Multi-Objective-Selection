function [fid, actual_path] = safe_fopen(preferred_path, mode, fallback_subdir)
%SAFE_FOPEN Try open a file; if fails, open in a writable temp fallback directory.
%
% [fid, actual_path] = safe_fopen(preferred_path, mode, fallback_subdir)
% - preferred_path: target path you want to write
% - mode: e.g. 'w'
% - fallback_subdir: e.g. 'morris_select_p_outputs'
%
% Returns:
% - fid: valid file id (never -1) or errors if both fail
% - actual_path: the actual path opened (preferred or fallback)

if nargin < 3 || isempty(fallback_subdir)
    fallback_subdir = 'morris_select_p_outputs';
end

% Ensure parent dir exists
pref_dir = fileparts(preferred_path);
if ~isempty(pref_dir) && ~exist(pref_dir, 'dir')
    mkdir(pref_dir);
end

% Try preferred path
fid = fopen(preferred_path, mode);
if fid ~= -1
    actual_path = preferred_path;
    return;
end

% Fallback to tempdir
fallback_dir = fullfile(tempdir, fallback_subdir);
if ~exist(fallback_dir, 'dir')
    mkdir(fallback_dir);
end
[~, name, ext] = fileparts(preferred_path);
fallback_path = fullfile(fallback_dir, [name ext]);

fid = fopen(fallback_path, mode);
if fid == -1
    error('SAFE_FOPEN:Failed', ...
        'Cannot open file for writing:\nPreferred: %s\nFallback:  %s\nClose programs locking files and check permissions.', ...
        preferred_path, fallback_path);
end

actual_path = fallback_path;
warning('Permission denied writing %s. Using fallback path: %s', preferred_path, fallback_path);
end
