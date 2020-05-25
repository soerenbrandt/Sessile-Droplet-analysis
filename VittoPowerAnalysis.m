function VittoPowerAnalysis
scale = 1/86.0006; % mm/pixel from reference image
init.vids = {'*.m4v','*.mov','*.avi','*.mp4'};

%---------------------------------------------------------------
%         Step 1 open video
%---------------------------------------------------------------
[name, path] = uigetfile({strjoin(init.vids,';'),'All Video Files'});
vid = [path, name];

v = VideoReader(vid);
videoLength = round(v.Duration*v.FrameRate);

%---------------------------------------------------------------
%         Step 2 select droplet to analyze and set the baseline
%---------------------------------------------------------------
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
viscircles(ax,center+rect(1:2), radius);
uiwait(h)

% Start new csv file
saveData(vid,radius*scale,CA(center,radius),'new');

%---------------------------------------------------------------
%         Step 3 finish video
%---------------------------------------------------------------
ImCount = 2;
w = waitbar(0,'Starting');
while hasFrame(v) && ImCount < 10000
    waitbar(ImCount/videoLength,w,['Frame ',num2str(ImCount),' of ',num2str(videoLength)]);
        
    currentImage = readFrame(v); % read the next video frame to analyze
    % find left edge
    [C1.center, C1.radius] = detectCircle(currentImage,[rect(1:2),rect(3)/2,rect(4)],setRadius);
    contactAngles = CA(C1.center,C1.radius); CAs.left = contactAngles.left;
    %find right edge
    [C2.center, C2.radius] = detectCircle(currentImage,[rect(1)+rect(3)/2,rect(2),rect(3)/2,rect(4)],setRadius);
    if size(C2.center,1) > 0; C2.center(1) = C2.center(1) + rect(3)/2; end
    contactAngles = CA(C2.center,C2.radius); CAs.right = contactAngles.right;
    radius = maxDistanceBetween(C1, C2);
    saveData(vid,radius*scale,CAs);
    
    ImCount = ImCount +1;
end
delete(w)

end

function dist = maxDistanceBetween(circle1, circle2)
try minX = circle1.center(1)-circle1.radius;
    maxX = circle2.center(1)+circle2.radius;
    dist = abs(maxX-minX)/2;
catch
    dist = NaN;
end
end

function [left, right, radius] = detectTwinCircle(Im,rect,estimate,baseline)
% Smoothen image for image processing
smoothIm = imopen(rgb2gray(imcrop(Im,rect)),strel('disk',15));
adjIm = imadjust(smoothIm);
% detect edge of droplet removing small isolated signals
dropletEdge = bwareaopen(edge(adjIm),10);
[Y, X] = find(dropletEdge); Y = size(adjIm,2)-Y; % coordinates of droplet edge
outlier = X*baseline(1)+baseline(2)>Y;% remove points below baseline
X(outlier) = []; Y(outlier) = [];
% separate approximate left and right edges
approxCenterX = size(adjIm,1)/2;
XL = X(X < approxCenterX); XR = X(X > approxCenterX);
YL = Y(X < approxCenterX); YR = Y(X > approxCenterX);
% fit circle path to droplet edges left and right
F = @(r, x0, y0 ,x)real(sqrt(r.^2 - (x-x0).^2))+y0; % circle function
c0 = [estimate size(adjIm,1)/2 size(adjIm,2)/2];
clb = [estimate*0.9 0 0]; cub = [estimate*1.1 size(adjIm,1) size(adjIm,2)];
FittedR = fit(XR,YR,F,'StartPoint',c0,'Lower',clb,'Upper',cub,'Robust','Bisquare');
FittedL = fit(XL,YL,F,'StartPoint',c0,'Lower',clb,'Upper',cub,'Robust','Bisquare');

% extract approximate radius
radius = (max(X) - min(X))/2;

% report circles
left.center = [FittedL.x0, size(adjIm,2)-FittedL.y0]; left.radius = FittedL.r;
right.center = [FittedR.x0, size(adjIm,2)-FittedR.y0]; right.radius = FittedR.r;
end

function [center, radius] = detectCircle(Im,rect,estimate)
% Smoothen image for image processing
smoothIm = imopen(rgb2gray(imcrop(Im,rect)),strel('disk',15));
adjIm = imadjust(smoothIm);
% find droplet edge and convert approximate drop by convex hull
dropAprox = bwconvhull(edge(adjIm));
%dropletRim = edge(imopen(rgb2gray(imcrop(Im,rect)),strel('disk',7))); % detect edge of droplet and create convex hull around points
% approximate background
th = 4;
bg = grayconnected(adjIm,1,1,th) + grayconnected(adjIm,size(adjIm,1),1,th) + grayconnected(adjIm,1,size(adjIm,2),th) + grayconnected(adjIm,size(adjIm,1),size(adjIm,2),th);
bg = bg > 0;
% highlight regions of image according to bg and dropAprox)
adjustments = ones(size(adjIm)); adjustments(dropAprox) = adjustments(dropAprox)+0.1; adjustments(bg>0) = adjustments(bg>0)-0.1;
finalAdjIm = double(adjIm).*adjustments/255;
% finalAdjIm = double(bwconvhull(dropletRim))*255;

[centers, radii, metric] = detectCircles(finalAdjIm,...
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
    CA.left = NaN; CA.right = NaN;
    if size(center) == [1, 2]; else; return; end 

    angle = @(m1,m2)rad2deg(rem(pi+atan((m1-m2)/(1+m1*m2)),pi));
    [x, y] =  linecirc(baseline(1),baseline(2),center(1),center(2), radius); % intersections circle with baseline
    for n = 1:length(x)
        tangentDirection = [x(n) y(n)] + [x(n)-center(1), y(n)-center(2)]*[0 -1; 1 0]; % tangent direction
        tangent = polyfit([x(n), tangentDirection(n)], [y(n), tangentDirection(n)], 1);
        
        if x(n) < center(1)
            CA.left = angle(tangent(n),baseline(1));
        else
            CA.right = angle(baseline(1),tangent(n));
        end
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