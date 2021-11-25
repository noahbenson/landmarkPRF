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
    u = linspace(-stimSize/2, stimSize/2, stimSize) * stimDegPerPx;
    [cols, rows] = meshgrid(u, u);
    tht = atan2(-rows, cols);
    ecc = hypot(cols, rows);
    ecclim = (ecc <= maxEccen);
    % We have a few categories of stimulus:
    % (1) Clockwise (30 s) and Counterclockwse (30 s)
    th0 = 0;
    thWidth = pi/2;
    nframesAng = 30*framesPerSecond;
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
    nframesEcc = 30*framesPerSecond;
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
    nframesBar = 20*framesPerSecond;
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
    % We also need a few pauses.
    blank15s = zeros([1,15*framesPerSecond], 'uint32');
    blank25s = zeros([1,25*framesPerSecond], 'uint32');
    % Now, make the images.
    seq = cat(2, blank15s, 1:size(masks, 3), blank25s);
    makeSeq(imgs, masks*255, seq, stimSize, framesPerSecond, outfile);
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
    % We want these eccentricities and these polar angles:
    wedgeWidths = [pi/16, pi/8, pi/4];
    wrMasks = zeros([stimSize, stimSize, 15], 'uint8');
    mx = wedgeWidths(1)/2;
    wrMasks(:,:,1) = (vrt < mx) & ecclim;
    wrMasks(:,:,2) = (hrz < mx) & ecclim;
    mx = wedgeWidths(2)/2;
    wrMasks(:,:,3) = (vrt < mx) & ecclim;
    wrMasks(:,:,4) = (hrz < mx) & ecclim;
    mx = wedgeWidths(3)/2;
    wrMasks(:,:,5) = (vrt < mx) & ecclim;
    wrMasks(:,:,6) = (hrz < mx) & ecclim;

    rings = {[0 0.5; 3.75 4.75; 7.5 8];
             [0 1;   3.5  4.5;  7   8];
             [0 2;   3    5;    6   8]};
    for ii = 1:3
        i0 = 3*(ii - 1) + 6;
        rr = rings{ii};
        wrMasks(:,:,1+i0) = (ecc >= rr(1,1)) & (ecc < rr(1,2));
        wrMasks(:,:,2+i0) = (ecc >= rr(2,1)) & (ecc < rr(2,2));
        wrMasks(:,:,3+i0) = (ecc >= rr(3,1)) & (ecc < rr(3,2));
    end
    % The wrMasks are:
    % 1: vertical meridian,   ±pi/8
    % 2: horizontal meridian, ±pi/8
    % 3: vertical meridian,   ±pi/4
    % 4: horizontal meridian, ±pi/4
    % 5: vertical meridian,   ±pi/2
    % 6: horizontal meridian, ±pi/2
    % 7: ring:  0-0.5°
    % 8: ring:  3.75°-4.25°
    % 9: ring:  7.5°-8°
    % 10: ring: 0°-1°
    % 11: ring: 3.5°-4.5°
    % 12: ring: 7°-8°
    % 13: ring: 0°-2°
    % 14: ring: 3°-5°
    % 15: ring: 6°-8°

    % Now make the stimulus!
    pauseFrames = zeros(1, 20*framesPerSecond, 'int32');
    ones4s = ones(1, 4*framesPerSecond, 'int32');
    zero4s = zeros(1, 4*framesPerSecond, 'int32');
    seq = cat(2, ...
        pauseFrames, ...
        7*ones4s, zero4s, 1*ones4s, zero4s, 8*ones4s, zero4s, 2*ones4s, zero4s, 9*ones4s, ...
        pauseFrames, ...
        10*ones4s, zero4s, 3*ones4s, zero4s, 11*ones4s, zero4s, 4*ones4s, zero4s, 12*ones4s, ...
        pauseFrames, ...
        13*ones4s, zero4s, 5*ones4s, zero4s, 14*ones4s, zero4s, 6*ones4s, zero4s, 15*ones4s, ...
        pauseFrames);
    makeSeq(imgs, wrMasks*255, seq, stimSize, framesPerSecond, outfile);
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