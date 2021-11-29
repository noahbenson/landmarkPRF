function lprf_runstim()
% Runs either the retinotopy experiment ("retinotopy") or the on-off wedge
% and ring experiment ("landmarks").

params = retCreateDefaultGUIParams;
params.experiment      = 'experiment from file';
params.fixation        = 'dot with grid';
params.modality        = 'fMRI';
params.stimSize        = 'max';
params.tr              = 1; % seconds
params.prescanDuration = 0; % seconds
params.calibration     = 'UWCHN_BOLDScreenLCD';
params.trigger         = 'scanner triggers computer';
params.triggerKey      = '5';
params.skipSyncTests   = true;
params.doEyelink       = false;
params.eyelinkIP       = '100.1.1.1';
params.quitProgKey     = 113; % 'q' key.
% If you want to save as a movie:
%params.movie = '/Users/nben/Desktop/retinotopy_stim.mov';

% These might not be necessary on all systems; for slower systems (i.e.,
% nben's MacBook Air), these are necessary to prevent Matlab from missing
% screen flips.
% waitBuffer: Pass control to PTB when the next flip is this far away (s).
params.waitBuffer = 0.1;
% waitPause: When there's at least twice the waitBuffer before the next
% screen flip, sleep this may seconds (then keep looping).
params.waitPause = 0.025;

% Instructions:
params.instructions = {
    'Thank you for participating in this vision experiment!';
    '';
    'During the experiment, please keep your eyes fixed on the';
    'green/red dot in the center of the screen.';
    'Whenever the dot changes color, press the response button.';
    '';
    'The experiment will start in a moment.';
    ''};
params.beginPrompt  = {
    'Thank you for participating in this vision experiment!';
    '';
    'During the experiment, please keep your eyes fixed on the';
    'green/red dot in the center of the screen.';
    'Whenever the dot changes color, press the response button.';
    '';
    'The experiment is about to begin.';
    'This text will disappear at that time.'};

fprintf('\n');
params.sub = input('Enter Subject ID: ', 's');
params.ses = input('Enter Session No: ', 's');

while true
    p = params;
    fprintf('\nStimulus Menu:\n');
    fprintf('  (1) Retinotopy\n');
    fprintf('  (2) Landmarks\n');
    fprintf('  (0) Quit\n');
    ch = input(' ? ', 's');
    if strcmp(ch, '1')
        expFile = 'stim_retinotopy';
        p.task = 'retinotopy';
    elseif strcmp(ch, '2')
        expFile = 'stim_landmark';
        p.task = 'landmarks';
    elseif strcmp(ch, '0')
        break;
    else
        fprintf('Bad input; try again.\n\n');
        continue;
    end

    % We load the experiment files from the directory containing this function.
    mfl = mfilename('fullpath');
    [mdir,~,~] = fileparts(mfl);
    p.loadMatrix = fullfile(mdir, [expFile '.mat']);
    
    % now set rest of the params
    p = setRetinotopyParams(p.experiment, p);
    % set response device
    p = setRetinotopyDevices(p);
    % Other params.
    p = retSetExperimentParams(p, p.experiment);
    p = retSetFixationParams(p, p.fixation);
    % go
    doRetinotopyScan(p);
end

