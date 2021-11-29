function lprf_genimgs(stimfile, outpath)
% Saves the sequence of stimulus images to the given output path, each with
% the filename 'frame%05d.png' at the appropriate resolution.
if ~exist(outpath, "dir")
    error("outpath doesn't exist: %s", outpath);
end
% Load the stimulus.
load(stimfile);
if ~exist('stimfile', 'var')
    error("stimfile does not contain 'stimulus' var");
end
% Walk through the seq/images.
for ii = 1:numel(stimulus.seq)
    outfl = fullfile(outpath, sprintf('frame%05d.png', ii-1));
    im = squeeze(stimulus.images(:,:,ii,:));
    im = draw_fixation(im);
    imwrite(im, outfl);
end

function im = draw_fixation(im)
    %% draw a polar grid to help the eye find fixation
    maxX = size(im,2);
    maxY = size(im,1);
    cenX = maxX ./ 2;
    cenY = maxY ./ 2;
    
	maxR = min([cenX, cenY]);
	rangeR = 0:100:maxR;
	rangeTheta = 0:30:330;
	
	for r = rangeR
        col = d.backColorIndex + 20;
		Screen('FrameOval', d.windowPtr, col, [cenX-r cenY-r cenX+r cenY+r]);
	end
		
	for th = rangeTheta
        col = d.backColorIndex + 20;
		[x y] = pol2cart(deg2rad(th), maxR);
		x = x + cenX;
		y = y + cenY;
		Screen('DrawLine', d.windowPtr, col, cenX, cenY, x, y);
	end
end
