function VittoPowerAnalysis
scale = 1/86.0006; % mm/pixel from reference image
init.vids = {'*.m4v','*.mov','*.avi','*.mp4'};

% open video
[name, path] = uigetfile({strjoin(init.vids,';'),'All Video Files'});
vid = [path, name];

v = VideoReader(vid);
videoLength = round(v.Duration*v.FrameRate);

%% select droplet to analyze and set the baseline
firstImage = readFrame(v);

h = figure; title(axes(h),'Select droplet'); % show first image and crop droplet
[Im, rect] = imcrop(firstImage);
delete(h)

h = figure; ax = axes(h); imshow(Im); hold(ax,'on'); title(ax,'Select left edge of drop') % Select baseline and size
c1 = ginput(1); title(ax,'Select right edge of drop') % Select left droplet edge
c2 = ginput(1); title(ax,'Select left baseline point') % Select right droplet edge
pos1 = ginput(1); title(ax,'Select right baseline point')% Select baseline point 1
pos2 = ginput(1); % Select baseline point 2
baseline = polyfit([pos1(1), pos2(1)], [pos1(2), pos2(2)], 1);
CA = @(c,r)contactAngle(c,r,baseline);
hold(ax,'off'); delete(h)

setRadius = abs(c1(1)-c2(1))/2;
[center, radius] = detectCircle(firstImage,rect,setRadius);
setRadius = min([setRadius,radius]); % update reference radius

h = figure; ax = axes(h); imshow(firstImage);
viscircles(ax,center+rect(1:2),radius);
uiwait(h)

% Start new file
saveData(vid,radius*scale,CA(center,radius),'new');

%% finish video
ImCount = 2;
w = waitbar(0,'Starting');
while hasFrame(v) && ImCount < 10000
    waitbar(ImCount/videoLength,w,['Frame ',num2str(ImCount),' of ',num2str(videoLength)]);
        
    currentImage = readFrame(v); % read the next video frame to analyze
    [center, radius] = detectCircle(currentImage,rect,setRadius);
    saveData(vid,radius*scale,CA(center,radius));
    
    ImCount = ImCount +1;
end
delete(w)

end



function [center, radius] = detectCircle(Im,rect,estimate)
% Calculate first image for verification
dropletRim = edge(imopen(rgb2gray(imcrop(Im,rect)),strel('disk',7))); % detect edge of droplet and create convex hull around points

[centers, radii, metric] = detectCircles(double(bwconvhull(dropletRim))*255,...
    [floor(0.9*estimate) ceil(1.1*estimate)],0.99); % find circle(s) that fit the convex hull

centers(metric<max(metric),:) = []; % remove poor quality fits
radii(metric<max(metric)) = [];

center = centers(radii == max(radii),:); % pick largest circle
radius = radii(radii == max(radii));
end

function saveData(vid,radius,CA, varargin)
[path,title] = fileparts(vid);
    
    fid=fopen([path,'/',title,'.csv'],'a');
    if fid < 0; error('The designated path is invalid'); end
    if strcmp('new',varargin)
        fprintf(fid,'%s, %s, %s\n','Radius /mm','CA left','CA right');
    end
    fprintf(fid,'%f, %f, %f\n',radius, CA.left, CA.right);
    fclose(fid);
end

function CA = contactAngle(center, radius, baseline)
    % validate input
    if size(center) == [1, 2]; else; CA.left = NaN; CA.right = NaN; return; end 

    angle = @(m1,m2)rad2deg(rem(pi+atan((m1-m2)/(1+m1*m2)),pi));
    [x, y] =  linecirc(baseline(1),baseline(2),center(1),center(2), radius); % intersections circle with baseline
    if length(x) == 2
        tangent1Direction = [x(1) y(1)] + [x(1)-center(1), y(1)-center(2)]*[0 -1; 1 0]; % tangent 1 direction
        tangent2Direction = [x(2) y(2)] + [x(2)-center(1), y(2)-center(2)]*[0 -1; 1 0]; % tangent 2 direction
        tangent1 = polyfit([x(1), tangent1Direction(1)], [y(1), tangent1Direction(2)], 1);
        tangent2 = polyfit([x(2), tangent2Direction(1)], [y(2), tangent2Direction(2)], 1);
        
        CA.left = angle(tangent1(1),baseline(1)); CA.right = angle(baseline(1),tangent2(1));
    elseif length(x) == 1
        tangent1Direction = [x(1) y(1)] + [x(1)-center(1), y(1)-center(2)]*[0 -1; 1 0]; % tangent 1 direction
        tangent1 = polyfit([x(1), tangent1Direction(1)], [y(1), tangent1Direction(2)], 1);
        
        if x(1) <= center(1)
            CA.left = angle(tangent1(1),baseline(1)); CA.right = NaN;
        else
            CA.left = NaN; CA.right = angle(baseline(1),tangent(1));
        end
    else
        CA.left = NaN; CA.right = NaN;
    end
    
    
end

function [centers, radii, metric] = detectCircles(Im,radius,varargin)
% finds circles in Im

% validate input
validateattributes(Im, {'numeric'},{'3d'})
validateattributes(radius, {'numeric'},{'nondecreasing','2d'})
if length(radius) == 1
	radius = [max([0, radius-15]) radius+15];
else
    validateattributes(radius, {'numeric'},{'size',[1,2]})
end

% check if Image is RGB and needs to be converted still
if size(Im,3) == 3
    Im = rgb2gray(Im);
end

% check input parameters
if strcmp('dark',varargin); polarity = 'dark'; else; polarity = 'bright'; end
if strcmp('TwoStage',varargin); method = 'TwoStage'; else; method = 'PhaseCode'; end
sensitivity = [varargin{cellfun(@(y)isa(y,'double'),varargin)}, 0.85];

% find circles
[centers, radii, metric] = imfindcircles(Im,radius,'ObjectPolarity',polarity,'Sensitivity',sensitivity(1), 'Method', method);

end