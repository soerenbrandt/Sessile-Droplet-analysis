function VittoPowerAnalysis
%scale = 1/86.0006; % mm/pixel from reference image
init.vids = {'*.m4v','*.mov','*.avi','*.mp4'};

%---------------------------------------------------------------
%         Step 1 open video
%---------------------------------------------------------------
[name, path] = uigetfile({strjoin(init.vids,';'),'All Video Files'});
vid = [path, name];

v = VideoReader(vid);
videoLength = round(v.Duration*v.FrameRate);

% request scale
response = inputdlg({'Number of pixels','per number of mm'},['Scale for ',name],1,{'86','1'});
scale = str2double(response{2})/str2double(response{1});

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
hold(ax,'off'); delete(h)

setRadius = abs(c1(1)-c2(1))/2;
drop = droplet(firstImage, rect, setRadius, baseline);
setRadius = min([setRadius,drop.Radius]); % update reference radius

h = plot(drop);
uiwait(h)

% Start new csv file
saveData(vid,drop.Radius*scale,drop.CA,drop.Height*scale, drop.Base*scale,'new');

%---------------------------------------------------------------
%         Step 3 finish video
%---------------------------------------------------------------
ImCount = 2;
w = waitbar(0,'Starting');
while hasFrame(v) && ImCount < 10000
    waitbar(ImCount/videoLength,w,['Frame ',num2str(ImCount),' of ',num2str(videoLength)]);
        
    currentImage = readFrame(v); % read the next video frame to analyze
    drop = droplet(currentImage, rect, setRadius, baseline);
    setRadius = min([setRadius,drop.Radius]); % update reference radius
    saveData(vid,drop.Radius*scale,drop.CA, drop.Height*scale, drop.Base*scale);
    
    ImCount = ImCount +1;
end
delete(w)

end


function saveData(vid,radius,CA, height, base, varargin)
[path,title] = fileparts(vid);
    
    fid=fopen([path,'/',title,'.csv'],'a');
    if fid < 0; error('The designated path is invalid'); end
    if strcmp('new',varargin)
        fprintf(fid,'%s, %s, %s, %s, %s\n','Radius /mm','CA left','CA right', 'Height /mm', 'Base radius /mm');
    end
    fprintf(fid,'%f, %f, %f, %f, %f\n',radius, CA.left, CA.right, height, base);
    fclose(fid);
end