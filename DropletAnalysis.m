function DropletAnalysisV5_1
FPS = 1;
warning('off')
init.vids = {'*.m4v','*.mov','*.avi','*.mp4'};

%---------------------------------------------------------------
%         Step 1 open video
%---------------------------------------------------------------
[name, path] = uigetfile({strjoin(init.vids,';'),'All Video Files'});
vid = [path, name];

v = VideoReader(vid);
videoLength = ceil(v.Duration);

% request scale
input = inputdlg({'Number of pixels','per number of mm','Show every # of fits'},['Scale for ',name],1,{'86','1','inf'});
scale = str2double(input{2})/str2double(input{1}); reportFreq = str2double(input{3});

%---------------------------------------------------------------
%         Step 2 select droplet to analyze and set the baseline
%---------------------------------------------------------------
firstImage = readFrame(v);

h = figure; title(axes(h),'Select droplet'); % show first image and crop droplet
[Im, rect] = imcrop(firstImage);
delete(h)

h = figure; ax = axes(h); imshow(Im,'InitialMagnification','fit'); 
hold(ax,'on'); title(ax,'Select left edge of drop') % Select baseline and size
c1 = ginput(1); title(ax,'Select right edge of drop') % Select left droplet edge
c2 = ginput(1); title(ax,'Select left baseline point') % Select right droplet edge
pos1 = ginput(1); title(ax,'Select right baseline point')% Select baseline point 1
pos2 = ginput(1); % Select baseline point 2
baseline = polyfit([pos1(1), pos2(1)], [pos1(2), pos2(2)], 1);
hold(ax,'off'); delete(h)

setRadius = abs(c1(1)-c2(1))/2;
drop = droplet(firstImage, rect, setRadius, baseline);
setRadius = min([setRadius,drop.Radius]); % update reference radius

h = plot(drop); % plot result of first frame and query

%---------------------------------------------------------------
%         Step 3 Query user desition to continue analysing (with or without
%                data video and continue
%---------------------------------------------------------------
response = questdlg('Continue analysis?','','Yes','Yes, with video','Cancel','Yes');

switch response
    case 'Yes'
    %---------------------------------------------------------------
    %         Step 3-1 finish video analysis with data output
    %---------------------------------------------------------------
    close(h)
    
    % Start new csv file
    saveData(vid, true, ...
                 'Time /s',0,'Radius /mm',drop.Radius*scale,'CA left', drop.CA.left, 'CA right', drop.CA.right, ...
                 'Height /mm', drop.Height*scale, 'Base radius /mm', drop.Base*scale, ...
                 'Volume /mm^3', drop.Volume*scale^3);

    % Finish video
    ImCount = 1;
    w = waitbar(0,'Starting');
    addTime = @(h,t)arrayfun(@(x)title(x,['Frame time: ',num2str(t),'s']),h.Children);
    timeRem = nan(1,100); estFrames = videoLength*FPS;
    while hasFrame(v) && ImCount < estFrames
        tic % start clock
        currentImage = readFrame(v); % read the next video frame to analyze

        if v.CurrentTime >= ImCount/FPS
        waitbar(ImCount/videoLength,w,['Second ',num2str(ImCount),' of ',num2str(videoLength),...
                                       ' (time rem: ',num2str(round((estFrames-ImCount)*nanmean(timeRem)/60)),'min)']);

        drop = droplet(currentImage, rect, setRadius, baseline);
        setRadius = min([setRadius,drop.Radius]); % update reference radius
        saveData(vid, false, ...
                 'Time /s',v.CurrentTime,'Radius /mm',drop.Radius*scale,'CA left', drop.CA.left, 'CA right', drop.CA.right, ...
                 'Height /mm', drop.Height*scale, 'Base radius /mm', drop.Base*scale, ...
                 'Volume /mm^3', drop.Volume*scale^3);

        if rem(ImCount,reportFreq) == 0
            h = plot(drop);
            addTime(h,v.CurrentTime);
        end

        ImCount = ImCount +1;
        timeRem(rem(ImCount,100)+1) = toc;
        end
    end
    delete(w)
    
    case 'Yes, with video'
    %---------------------------------------------------------------
    %         Step 3-2 finish video analysis including output video
    %---------------------------------------------------------------
    % Preallocate movie
    editedFrames(ceil(v.Duration)) = getframe(h);
    editedFrames = fliplr(editedFrames);
    close(h)
    
    % Start new csv file
    saveData(vid, true, ...
                 'Time /s',0,'Radius /mm',drop.Radius*scale,'CA left', drop.CA.left, 'CA right', drop.CA.right, ...
                 'Height /mm', drop.Height*scale, 'Base radius /mm', drop.Base*scale, ...
                 'Volume /mm^3', drop.Volume*scale^3);
    
    % Start new csv file
    saveData(vid, true, ...
                 'Time /s',0,'Radius /mm',drop.Radius*scale,'CA left', drop.CA.left, 'CA right', drop.CA.right, ...
                 'Height /mm', drop.Height*scale, 'Base radius /mm', drop.Base*scale, ...
                 'Volume /mm^3', drop.Volume*scale^3);

    % Finish video
    ImCount = 1;
    w = waitbar(0,'Starting');
    addTime = @(h,t)arrayfun(@(x)title(x,['Frame time: ',num2str(t),'s']),h.Children);
    timeRem = nan(1,100); estFrames = videoLength*FPS;
    while hasFrame(v) && ImCount < videoLength*FPS
        tic % start clock
        currentImage = readFrame(v); % read the next video frame to analyze

        if v.CurrentTime >= ImCount/FPS
        waitbar(ImCount/videoLength,w,['Second ',num2str(ImCount),' of ',num2str(videoLength),...
                                       ' (time rem: ',num2str(round((estFrames-ImCount)*nanmean(timeRem)/60)),'min)']);

        drop = droplet(currentImage, rect, setRadius, baseline);
        setRadius = min([setRadius,drop.Radius]); % update reference radius
        saveData(vid, false, ...
                 'Time /s',v.CurrentTime,'Radius /mm',drop.Radius*scale,'CA left', drop.CA.left, 'CA right', drop.CA.right, ...
                 'Height /mm', drop.Height*scale, 'Base radius /mm', drop.Base*scale, ...
                 'Volume /mm^3', drop.Volume*scale^3);

        if rem(ImCount,reportFreq) == 0
            h = plot(drop);
            addTime(h,v.CurrentTime);
        end
        
        % Plot movie frame
        movFrame = plot(drop);
        set(movFrame, 'Visible', 'off');
        
        addTime(movFrame,v.CurrentTime);
%         addDataToFrame(movFrame, ...
%                  'Radius /mm',drop.Radius*scale,'CA left', drop.CA.left, 'CA right', drop.CA.right, ...
%                  'Height /mm', drop.Height*scale, 'Base radius /mm', drop.Base*scale, ...
%                  'Volume /mm^3', drop.Volume*scale^3);
        editedFrames(ImCount+1) = getframe(movFrame);
        
        ImCount = ImCount +1;
        timeRem(rem(ImCount,100)+1) = toc;
        end
    end
    
    waitbar(1,w,strcat('writing video file'));
    % Prepare new videoFile for writing
    [path,name] = fileparts(vid);
    newVideo = VideoWriter([path,'/',name,' - analyzed'],'MPEG-4');%'Uncompressed AVI');
    newVideo.FrameRate = FPS;

    % Clean out empty editedFrames
    editedFrames(arrayfun(@(x)isempty(x.cdata),editedFrames)) = [];

    % Write video file
    open(newVideo)
    writeVideo(newVideo,editedFrames);
    close(newVideo)

    % delete last waitbar
    delete(w)

    otherwise % do nothing
        close(h)
        return
end

end


% function addDataToFrame(h,varargin)
% arrayfun(@(n)validateattributes(varargin{n},{'char'},{}),1:2:numel(varargin))
% arrayfun(@(n)validateattributes(varargin{n},{'numeric'},{}),2:2:numel(varargin))
% if isempty(varargin)
%     error('Expects at least one string-value input pair.')
% elseif rem(numel(varargin),2)
%     error('Variable input expected to be string-value pairs.')
% end
% try varargin{arrayfun(@(n)isa(varargin{n},'numeric') && isnan(varargin{n}),1:numel(varargin))} = [];
% catch
% end
% 
% dim = [0.5 0.5 0.3 0.3];
% str = arrayfun(@(n)[varargin{n},' = ',num2str(varargin{n+1})],1:2:numel(varargin),'uni',0);
% annotation(h,'textbox',dim,'String',str,'FitBoxToText','on','FontSize',10,'Interpreter','latex');
% 
% end

function saveData(vid, new, varargin)
arrayfun(@(n)validateattributes(varargin{n},{'char'},{}),1:2:numel(varargin))
arrayfun(@(n)validateattributes(varargin{n},{'numeric'},{}),2:2:numel(varargin))
if isempty(varargin)
    error('Expects at least one string-value input pair.')
elseif rem(numel(varargin),2)
    error('Variable input expected to be string-value pairs.')
end
try varargin{arrayfun(@(n)isa(varargin{n},'numeric') && isnan(varargin{n}),1:numel(varargin))} = [];
catch
end
try varargin(cellfun(@(x)~isa(x,'char') && any(isnan(x)), varargin)) = {''}; % replace NaN with white space
catch
    disp('now')    
end

[path,title] = fileparts(vid);
    
    fid=fopen([path,'/',title,'.csv'],'a');
    if fid < 0; error('The designated path is invalid'); end
    if new
        fprintf(fid,['%s',repmat(', %s',1,numel(varargin)/2-1),'\n'],varargin{1:2:end});
    end
    fprintf(fid,['%f',repmat(',%f',1,numel(varargin)/2-1),'\n'],varargin{2:2:end});
    fclose(fid);
end
