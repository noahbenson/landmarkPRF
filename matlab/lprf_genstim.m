function genStimulus(varargin)
%% Generate the Retinotopy Stimulus for VistaDisp using the HCP stimulus.
%  This script generates the stimulus, as understood by the vistadisp
%  repository (https://github.com/noahbenson/vistadisp), for retinotopic
%  mapping and retinotopic landmark stimulation. The script downloads the
%  HCP stimulus data from kendrickkay.net if needed, and writes all
%  stimulus files to the directory containing the script.

    %% Some constants for this script.
    %  Screen data. These data are regarding the Screen for the UW CHN's
    %  Cambridge Research Systems BOLDScreen 32" LCD; they were measured by
    %  the UW CHN in 2021.
    stimDegPerPx    = 0.0157;
    screenSize      = 1080;
    %  Stimulus Details. All stimulus files will have the following
    %  frame-rate; additionally, the max eccentricity of the stimulus can
    %  be used to reduce the image size slightly.
    framesPerSecond = 10;
    maxEccen        = 8;
    %  We set the stimSize to be based on the max eccentricity.
    stimSize = min(screenSize, ceil(maxEccen * 2 / stimDegPerPx));
    
    % Check the input is either 0 or 1 item.
    if numel(varargin) > 1, error("only 1 argument may be given"); end
    patterns = varargin;

    mfl = mfilename('fullpath');
    [mdir,~,~] = fileparts(mfl);
    outfile = fullfile(mdir, 'stim_retinotopy.mat');
    if ~exist(outfile, 'file')
        patterns = {loadPatterns(mdir, patterns{:})};
        genRetinotopyStimulus(patterns{1}, ...
                              stimSize, ...
                              stimDegPerPx, ...
                              framesPerSecond, ...
                              maxEccen, ...
                              outfile);
    end
    
    % For the landmark stimulus, we can just make our own masks.
    outfile = fullfile(mdir, 'stim_landmark.mat');
    if ~exist(outfile, 'file')
        patterns = {loadPatterns(mdir, patterns{:})};
        genLandmarkStimulus(patterns{1}, ...
                            stimSize, ...
                            stimDegPerPx, ...
                            framesPerSecond, ...
                            maxEccen, ...
                            outfile);
    end
end

% This function makes the stimulus sequence given a set of images to
% flicker in the masks, a set of masks, a mask sequence, the stimulus
% size, frames per second, and the output file.
function makeSeq(imgs, msks, seq0, stimSize, framesPerSecond, outfile)
    fprintf('Making stimulus file %s...\n', outfile);
    % We can make the seq now; it's the same for both retinotopy runs.
    nframes = numel(seq0);
    seqtiming = ((1:nframes) - 1) / framesPerSecond;
    % The output frames we will write into.
    frames = 128 * ones([stimSize stimSize nframes 3], 'uint8');
    fprintf('  - Collecting frames...\n');
    scale = stimSize / size(imgs,1);
    nimgs = size(imgs,4);
    for ii = 1:nframes
        s = seq0(ii);
        if s == 0
            frames(:,:,ii,:) = 128;
        else
            imno = mod(ii - 1, nimgs) + 1; %randi(nimgs);
            hcpim = squeeze(imgs(:,:,:,imno));
            if size(hcpim,1) ~= stimSize
                hcpim = imresize(hcpim, scale);
            end
            mask = msks(:,:,s);
            if stimSize ~= size(mask,1)
                mask = imresize(mask, scale);
            end
            mask = uint8(round(mask / 255));
            im = uint8(128 .* (1 - mask)) + uint8(hcpim .* mask);
            frames(:,:,ii,:) = im;
        end
    end
    % Save it.
    stimulus = [];
    stimulus.images = frames;
    stimulus.seq = (1:nframes)';
    stimulus.seqtiming = seqtiming(:)';
    stimulus.cmap = [0:255; 0:255; 0:255]' / 255;
    fprintf('  - Saving file...\n');
    save(outfile, 'stimulus', '-v7.3');
end

% Generate the Retinotopy Stimulus sequence.
function genRetinotopyStimulus(imgs, stimSize, stimDegPerPx, ...
                               framesPerSecond, maxEccen, outfile)
    fprintf('Generating Retinotopy Masks...\n');
    masks = genRetinotopyMasks(stimSize, stimDegPerPx, ...
                               maxEccen, framesPerSecond, 30, 20);
    % We also need a few pauses.
    blank15s = zeros([1,15*framesPerSecond], 'uint32');
    blank25s = zeros([1,25*framesPerSecond], 'uint32');
    % Now, make the images.
    seq = cat(2, blank15s, 1:size(masks, 3), blank25s);
    makeSeq(imgs, masks*255, seq, stimSize, framesPerSecond, outfile);
end
function masks = genRetinotopyMasks(stimSize, stimDegPerPx, ...
                                    maxEccen, framesPerSecond, ...
                                    secPerWR, secPerBar)
    u = linspace(-stimSize/2, stimSize/2, stimSize) * stimDegPerPx;
    [cols, rows] = meshgrid(u, u);
    tht = atan2(-rows, cols);
    ecc = hypot(cols, rows);
    ecclim = (ecc <= maxEccen);
    % We have a few categories of stimulus:
    % (1) Clockwise (30 s) and Counterclockwse (30 s)
    th0 = 0;
    thWidth = pi/2;
    nframesAng = secPerWR*framesPerSecond;
    masksCCW = zeros([stimSize, stimSize, nframesAng], 'uint8');
    for ii = 1:nframesAng
        % Clockwise:
        th = th0 + 2*pi * (ii-1)/nframesAng;
        masksCCW(:,:,ii) = genWedgeMask(tht, ecclim, th, thWidth/2);
    end
    % Counterclockwise:
    masksCW = flip(masksCCW, 2);
    % (2) Expanding (30 s) and Contracting (30 s).
    %     We grow the ring to half / contract it to half then start moving
    %     the inner bound out / outer bound in.
    nframesEcc = secPerWR*framesPerSecond;
    masksEXP = zeros([stimSize, stimSize, nframesEcc], 'uint8');
    for ii = 1:nframesEcc
        frac = (ii-1) / (nframesEcc - 1);
        if frac < 1/3
            ecc1 = maxEccen/2 * (frac * 3);
            ecc0 = 0;
        elseif frac < 2/3
            ecc0 = (frac - 1/3)*1.5 * maxEccen;
            ecc1 = ((frac - 1/3)*1.5 + 0.5) * maxEccen;
        else
            ecc0 = (frac - 1/3)*1.5 * maxEccen;
            ecc1 = maxEccen;
        end
        masksEXP(:,:,ii) = genRingMask(ecc, ecc0, ecc1);
    end
    masksCON = flip(masksEXP, 3);
    % (3) Left-Right and Right-Left bars; we use a similar strategy as the
    %     eccentricity mapping, but we will fill only 1/4 of the screen.
    nframesBar = secPerBar*framesPerSecond;
    masksLR = zeros([stimSize, stimSize, nframesBar], 'uint8');
    for ii = 1:nframesBar
        x1 = maxEccen * (2.5 * ii/nframesBar - 1);
        x0 = maxEccen * (2.5 * ii/nframesBar - 1.5);
        masksLR(:,:,ii) = genBarMask(cols, ecclim, x0, x1);
    end
    masksRL = flip(masksLR, 2);
    % (4) Down-Up and Up-Down bars.
    masksUD = permute(masksLR, [2,1,3]);
    masksDU = flip(masksUD, 1);
    % That's all the masks; we'll be going through each of them once with
    % some breaks in between; go ahead and setup the sub-sequences and the
    % stack of masks for this.
    masks = cat(3, ...
                masksCON, masksRL, ...
                masksCCW, masksDU, ...
                masksLR,  masksCW, ...
                masksUD,  masksEXP);
end
function mask = genWedgeMask(theta_im, ecclim, theta, width)
    mask = ecclim & (abs(mod(theta_im - theta - pi, 2*pi) - pi) < width);
end
function mask = genRingMask(ecc_im, ecc0, ecc1)
    mask = (ecc_im >= ecc0) & (ecc_im < ecc1);
end
function mask = genBarMask(x_im, ecclim, x0, x1)
    mask = ecclim & (x_im >= x0) & (x_im < x1);
end

% Generate the Landmark Stimulus sequence.
function genLandmarkStimulus(imgs, stimSize, stimDegPerPx, ...
                             framesPerSecond, maxEccen, outfile)
    fprintf('Generating Landmark Masks...\n');
    u = linspace(-stimSize/2, stimSize/2, stimSize) * stimDegPerPx;
    [cols, rows] = meshgrid(u, u);
    tht = atan2(-rows, cols);
    hrz = min(abs(tht), fliplr(abs(tht)));
    vrt = hrz';
    ecc = hypot(cols, rows);
    ecclim = (ecc <= maxEccen);
    % We need to generate three experiments:
    % (1) fast retinotopy
    % (2) landmarks
    % (3) subdivisions
    % (4, 5, 6) repeat the above
    % Each is separated by 20 seconds with 15 seconds at the start
    % and 25 seconds at the end of the experiment (2:20 of blank).
    % Each experiment lasts 30 seconds (3:00 of experiment).
    % Start with retinotopy: we use the above function...
    ret_masks = genRetinotopyMasks(stimSize, stimDegPerPx, ...
                                   maxEccen, framesPerSecond, 4.5, 3);
    ret_seq = 1:size(ret_masks,3);
    % Landmarks next;
    % The wrMasks are:
    % 1: vertical meridian,   ±pi/16
    % 2: horizontal meridian, ±pi/16
    % 3: ring:  0-0.5°
    % 4: ring:  3.75°-4.25°
    % 5: ring:  7.5°-8°
    wedgeWidth = pi/16;
    wr_masks = zeros([stimSize, stimSize, 5], 'uint8');
    wr_masks(:,:,1) = (vrt < wedgeWidth/2) & ecclim;
    wr_masks(:,:,2) = (hrz < wedgeWidth/2) & ecclim;
    rr = [0 0.5; 3.75 4.75; 7.5 8];
    wr_masks(:,:,3) = (ecc >= rr(1,1)) & (ecc < rr(1,2));
    wr_masks(:,:,4) = (ecc >= rr(2,1)) & (ecc < rr(2,2));
    wr_masks(:,:,5) = (ecc >= rr(3,1)) & (ecc < rr(3,2));
    seq_long = ones([1, 4*framesPerSecond], 'int32');
    seq_shrt = zeros([1, 2.5*framesPerSecond], 'int32');
    wr_seq = cat(2, ...
                 1*seq_long, seq_shrt, 2*seq_long, seq_shrt, ...
                 3*seq_long, seq_shrt, 4*seq_long, seq_shrt, ...
                 5*seq_long);
    % Now, make the subdivisions:
    % They are 2 sec on, 2 sec off.
    sub_masks = zeros([stimSize, stimSize, 8], 'uint8');
    sub_masks(:,:,1) = (mod(tht, pi) <  pi/2) & ecclim;
    sub_masks(:,:,2) = (mod(tht, pi) >= pi/2) & ecclim;
    sub_masks(:,:,3) = (mod(tht, pi/2) <  pi/4) & ecclim;
    sub_masks(:,:,4) = (mod(tht, pi/2) >= pi/4) & ecclim;
    sub_masks(:,:,5) = (ecc <  maxEccen/2);
    sub_masks(:,:,6) = (ecc >= maxEccen/2) & ecclim;
    sub_masks(:,:,7) = (mod(ecc, maxEccen/2) <  maxEccen/4) & ecclim;
    sub_masks(:,:,8) = (mod(ecc, maxEccen/2) >= maxEccen/4) & ecclim;
    seq_2s = ones([1, 2*framesPerSecond], 'int32');
    sub_seq = cat(2, ...
                  1*seq_2s, 0*seq_2s, 2*seq_2s, 0*seq_2s, ...
                  3*seq_2s, 0*seq_2s, 4*seq_2s, 0*seq_2s, ...
                  5*seq_2s, 0*seq_2s, 6*seq_2s, 0*seq_2s, ...
                  7*seq_2s, 0*seq_2s, 8*seq_2s, 0*seq_2s);
    % Now make the stimulus!
    ii = wr_seq > 0;
    wr_seq(ii) = wr_seq(ii) + size(ret_masks, 3);
    ii = sub_seq > 0;
    sub_seq(ii) = sub_seq(ii) + size(ret_masks, 3) + size(wr_masks, 3);
    seq_15s = zeros([1, 15*framesPerSecond], 'int32');
    seq_25s = zeros([1, 15*framesPerSecond], 'int32');
    seq_20s = zeros([1, 15*framesPerSecond], 'int32');
    seq = [seq_15s ...
           ret_seq seq_20s wr_seq seq_20s sub_seq seq_20s ...
           ret_seq seq_20s wr_seq seq_20s sub_seq ...
           seq_25s];
    masks = cat(3, ret_masks, wr_masks, sub_masks);
    makeSeq(imgs, masks*255, seq, stimSize, framesPerSecond, outfile);
end

% Load the HCP Stimulus
function patterns = loadPatterns(mdir, varargin)
    if numel(varargin) == 1
        patterns = varargin{1};
        return;
    end
    hcpstimfl = fullfile(mdir, 'HCP_stimdata.mat');
    if ~exist(hcpstimfl, 'file')
        fprintf('Downloading HCP stimulus file from kendrickkay.net...\n')
        url = 'http://kendrickkay.net/analyzePRF/stimuli.mat';
        dat = webread(url);
        fl = fopen(hcpstimfl, 'w');
        fwrite(fl, dat);
        fclose(fl);
    end
    fprintf('Loading HCP pattern images...\n');
    load(hcpstimfl, 'patterns');
end