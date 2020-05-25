classdef droplet
    % Comes with VittoPowerAnalysisV2 to interpret droplet shape and size
    properties
        Shape % array of ellipses that describe the shape of the drop
        CA % contact angles with the baseline of the droplet (left and right)
    end
    properties (Dependent)
        Radius % radius dimension of the droplet in px
        Base % radius dimension of the droplet in px
        Height % height dimension of the droplet in px
        Volume % approximate volume of the droplet
    end
    properties (SetAccess = immutable)
        estimate % radius estimate at initialization
        baseline % defined baseline in the image
        Image % original image given
    end
    
    methods
        function obj = droplet(Im, rect, estimate, baseline) % calculates properties of the droplet including size and CA
            %---------------------------------------------------------------
            %         Step 1 Initialize immutable properties
            %---------------------------------------------------------------
            
            obj.Image = imcrop(Im,rect);
            obj.baseline = baseline;
            obj.estimate = estimate;
            
            quality = gaugeQuality(obj);
            
            %---------------------------------------------------------------
            %         Step 2 Aproximate shape of the droplet based on image
            %         quality
            %---------------------------------------------------------------
            switch quality
                case 1
                    [obj.Shape, obj.CA] = doubleEllipseAprox(obj);
                case 2
                    [obj.Shape, obj.CA] = sphericalAprox(obj);
                otherwise
                    [obj.Shape, obj.CA] = sphericalContactAprox(obj);
            end
        end
        function fig = plot(obj)
            %---------------------------------------------------------------
            %         Step 0 Define regions for each part of the shape
            %---------------------------------------------------------------
            angle = @(m1,m2)pi-rem(pi+atan((m1-m2)/(1+m1*m2)),pi);
            left = (-1.35*pi : 0.01 : -0.55*pi) + angle(0,obj.baseline(1));
            middle = (-0.55*pi : 0.01 : -0.45*pi) + angle(0,obj.baseline(1));
            right = (-0.45*pi : 0.01 : 0.35*pi) + angle(0,obj.baseline(1));
            
            %---------------------------------------------------------------
            %         Step 1 Create ellipses
            %---------------------------------------------------------------
            R = @(phi)[ cosd(phi) sind(phi); -sind(phi) cosd(phi) ]; % rotation matrix in 2D
            X = @(a,X0,range)  a * cos(range) + X0; Y = @(b,Y0,range) b * sin(range) + Y0; % functions plotting X,Y-values of ellipse
            rotEll = @(ell,range) R(ell.phi) * [X(ell.a,ell.X0,range+deg2rad(ell.phi)); Y(ell.b,ell.Y0,range+deg2rad(ell.phi))];
            E1 = rotEll(obj.Shape(1),left);
            Em = rotEll(obj.Shape(2),middle);
            E2 = rotEll(obj.Shape(3),right);
            
            %---------------------------------------------------------------
            %         Step 2 Plot
            %---------------------------------------------------------------
            fig = figure; h = axes(fig);
            imshow(obj.Image,'InitialMagnification','fit')
            hold(h,'on')
                plot(h,E1(1,:),E1(2,:),'r-'); % plot shape - left
                plot(h,Em(1,:),Em(2,:),'r-'); % plot shape - middle
                plot(h,E2(1,:),E2(2,:),'b-'); % plot shape - right
                
                lr = [0.05*size(obj.Image,2) 0.95*size(obj.Image,2)];
                plot(h,linspace(lr(1),lr(2),2),obj.baseline(1)*linspace(lr(1),lr(2),2)+obj.baseline(2),'w-'); % plot baseline
                plot(h,linspace(lr(1),lr(2),2),obj.baseline(1)*linspace(lr(1),lr(2),2)+obj.baseline(2)-obj.Height,'w--'); % plot height
                
                % plot contact angles
                L1 = obj.CA.leftPin - 25*[obj.baseline(2) 1]*R(-obj.CA.left)/norm([obj.baseline(2) 1]); 
                L2 = L1 + 75*[obj.baseline(2) 1]*R(-obj.CA.left)/norm([obj.baseline(2) 1]);
                line(h,[L1(1) L2(1)],[L1(2) L2(2)],'Color','w','LineStyle','-');
                R1 = obj.CA.rightPin - 25*[obj.baseline(2) 1]*R(180+obj.CA.right)/norm([obj.baseline(2) 1]); 
                R2 = R1 + 75*[obj.baseline(2) 1]*R(180+obj.CA.right)/norm([obj.baseline(2) 1]);
                line(h,[R1(1) R2(1)],[R1(2) R2(2)],'Color','w','LineStyle','-');
                
                % plot diameter
                [Pl, Pr] = obj.pointsFurthestApart;
                plot([Pl(1) Pr(1)], [Pl(2) Pr(2)], 'Color', 'w')
            hold(h,'off')
        end % CHANGE
        function [shape, CA] = doubleEllipseAprox(obj,varargin)
            % Separately fits ellipses to left and right edge
            
            %---------------------------------------------------------------
            %         Step 1 Detect edge of the droplet
            %---------------------------------------------------------------
            H = outlineDrop(obj); % conversion to grayscale, edge detection, fills in small gaps in edge, removes small false positives, adds bottom line
            I = imfill(H,'holes'); I = bwareaopen(I,round(sum(I(:))/2)); % fill enclosed area and remove remaining small outliers
                I = bwmorph(I,'remove'); 
            line = @(x)round(obj.baseline(2) + obj.baseline(1)*x);    
            I(sub2ind(size(H),line(1:size(H,1)), 1:size(H,1))) = 0; % I(size(I,1),:) = 0;
                
            %---------------------------------------------------------------
            %         Step 2 Separate left and right side of the edge
            %---------------------------------------------------------------
            [Y,X] = find(I); leftMax = min(X(Y == min(Y))); rightMin = max(X(Y == min(Y)));
            left = struct('X',X(X < leftMax),'Y',Y(X < leftMax));
            right = struct('X',X(X > rightMin),'Y',Y(X > rightMin));
            
            % remove reflection underneath baseline
            distToBase = @(x,y)abs(obj.baseline(1).*x-y+obj.baseline(2))./sqrt(obj.baseline(1)^2+1);
            
            [~,leftBottom] = min(distToBase(left.X,left.Y)); 
                left.X = left.X(left.Y<left.Y(leftBottom));
                left.Y = left.Y(left.Y<left.Y(leftBottom));
            [~,rightBottom] = min(distToBase(right.X,right.Y));
                right.X = right.X(right.Y<right.Y(rightBottom));
                right.Y = right.Y(right.Y<right.Y(rightBottom));
                
            %---------------------------------------------------------------
            %         Step 3 Fit ellipses to left and right side
            %---------------------------------------------------------------    
            ellL = fit_ellipse(left.X,left.Y);
            ellR = fit_ellipse(right.X,right.Y);
            
            %---------------------------------------------------------------
            %         Step 4 Report droplet shape and contact angles
            %---------------------------------------------------------------
            E1 = droplet.createEllipse([ellL.a, ellL.b], [ellL.X0, ellL.Y0], ellL.phi);
            E2 = droplet.createEllipse([ellR.a, ellR.b], [ellR.X0, ellR.Y0], ellR.phi);
            
            % find middle
            circle.radius = max(distToBase(X,Y));
            m = floor(median(find(distToBase(X,Y) == circle.radius)));
            circle.center = [X(m),Y(m)] + circle.radius*([1 obj.baseline(1)]*[0 1; -1 0]/norm([obj.baseline(1) 1]));
            Em = droplet.createEllipse(circle.radius, circle.center, 0);
            
            shape = [E1, Em, E2];
            
            leftCA = contactAngle(obj,E1,-1); rightCA = contactAngle(obj,E2,1);
            CA = struct('left',leftCA.CA,'leftPin',leftCA.pin,...
                         'right',rightCA.CA,'rightPin',rightCA.pin);
        end
        function [shape, CA] = sphericalAprox(obj)
            % Use Hough transform to detect circle and keep
            % largest, best fit
            
            %---------------------------------------------------------------
            %         Step 1 Smoothen image for image processing
            %---------------------------------------------------------------
            smoothIm = imopen(rgb2gray(obj.Image),strel('disk',15));
            adjIm = imsharpen(imadjust(smoothIm));
            
            %---------------------------------------------------------------
            %         Step 1 Find best circle fit in whole image
            %---------------------------------------------------------------
            [centers, radii, metric] = imfindcircles(adjIm,[floor(0.9*obj.estimate) ceil(1.1*obj.estimate)],'Sensitivity',0.99,'ObjectPolarity','dark');
            
            centers(metric<max(metric),:) = []; % remove poor quality fits
            radii(metric<max(metric)) = [];

            circle.center = mean(centers(radii == max(radii),:),1); % pick largest circle
            circle.radius = mean(radii(radii == max(radii)));
            
            %---------------------------------------------------------------
            %         Step 2 Report droplet shape and contact angles
            %---------------------------------------------------------------
            E1 = droplet.createEllipse(circle.radius, circle.center, 0);
            shape = [E1, E1, E1];
            
            leftCA = contactAngle(obj,E1,-1); rightCA = contactAngle(obj,E1,1);
            CA = struct('left',leftCA.CA,'leftPin',leftCA.pin,...
                         'right',rightCA.CA,'rightPin',rightCA.pin);
        end
        function [shape, CA] = sphericalContactAprox(obj)
            % Separately image the right half as dark using Hough transform
            % to get a better approximation
            
            %---------------------------------------------------------------
            %         Step 0 Smoothen image for image processing
            %---------------------------------------------------------------
            smoothIm = imopen(rgb2gray(obj.Image),strel('disk',15));
            adjIm = imsharpen(imadjust(smoothIm));
            
            %---------------------------------------------------------------
            %         Step 1 Find left edge using whole image
            %---------------------------------------------------------------
            [centers, radii, metric] = droplet.detectCircles(adjIm,[floor(0.9*obj.estimate) ceil(1.1*obj.estimate)],0.99);
            
            centers(metric<max(metric),:) = []; % remove poor quality fits
            radii(metric<max(metric)) = [];

            circle.center = centers(radii == max(radii),:); % pick largest circle
            circle.radius = radii(radii == max(radii));
            E1 = droplet.createEllipse(circle.radius, circle.center, 0);
            
            %---------------------------------------------------------------
            %         Step 2 Find right edge using right half of image
            %---------------------------------------------------------------
            imDim = size(adjIm); rightHalf = [ceil(imDim(1)/2) 1 imDim(1)-ceil(imDim(1)/2) imDim(2)];
            [centers, radii, metric] = droplet.detectCircles(imcrop(adjIm,rightHalf),[floor(0.9*obj.estimate) ceil(1.1*obj.estimate)],0.99,'dark');
            if size(centers,1) > 0; centers(:,1) = centers(:,1) + ceil(imDim(1)/2); end
            
            centers(metric<max(metric),:) = []; % remove poor quality fits
            radii(metric<max(metric)) = [];

            circle.center = centers(radii == max(radii),:); % pick largest circle
            circle.radius = radii(radii == max(radii));
            E2 = droplet.createEllipse(circle.radius, circle.center, 0);
            
            %---------------------------------------------------------------
            %         Step 3 Report droplet shape and contact angles
            %---------------------------------------------------------------
            shape = [E1, E1, E2];
            
            leftCA = contactAngle(obj,E1,-1); rightCA = contactAngle(obj,E2,1);
            CA = struct('left',leftCA.CA,'leftPin',leftCA.pin,...
                         'right',rightCA.CA,'rightPin',rightCA.pin);
        end
        function radius = get.Radius(obj) % estimates the radius of the droplet
            [Pl, Pr] = obj.pointsFurthestApart;
            
            radius = sqrt((Pr(1) - Pl(1))^2 + (Pr(2) - Pl(2))^2)/2;
        end
        function base = get.Base(obj) % calculates the radius of the base of the droplet
            try 
                left = obj.CA.leftPin;
                right = obj.CA.rightPin;
                base = sqrt((right(1)-left(1))^2+(right(2)-left(2))^2)/2;
            catch
                base = []; return
            end
        end
        function height = get.Height(obj) % calculates the height of the droplet over the baseline
            try shapeAtHeight = obj.Shape(2); % Middle circle describes the height best
                distToBase = abs(obj.baseline(1).*shapeAtHeight.X0-shapeAtHeight.Y0+obj.baseline(2))./sqrt(obj.baseline(1)^2+1);
                height = distToBase + shapeAtHeight.a;
            catch
                height = [];
            end
        end
        function volume = get.Volume(obj) % approximates the volume of the droplet over the baseline
            volume = nan;
            if gaugeQuality(obj) ~= 1
                return
            end
            
            %---------------------------------------------------------------
            %         Step 1 Find center point on baseline and ensure it is
            %                viable
            %---------------------------------------------------------------
            center = [(obj.CA.leftPin(1) + obj.CA.rightPin(1))/2 ...
                          (obj.CA.leftPin(2) + obj.CA.rightPin(2))/2];
            if any(isnan(center))
                return
            end
            
            %---------------------------------------------------------------
            %         Step 2 Outline shape of the droplet
            %---------------------------------------------------------------
            H = outlineDrop(obj);
            I = imfill(H,'holes'); I = bwareaopen(I,round(sum(I(:))/2)); % fill enclosed area and remove remaining small outliers
                I = bwmorph(I,'remove'); I(size(I,1),:) = 0;
            [Y,X] = find(I);
            
            %---------------------------------------------------------------
            %         Step 3 Rotate center point and outline
            %---------------------------------------------------------------
            R = @(phi)[ cos(phi) sin(phi); -sin(phi) cos(phi) ]; % rotation matrix in 2D
            angle = @(m1,m2)pi-rem(pi+atan((m1-m2)/(1+m1*m2)),pi);
            
            centerp = center * R(angle(obj.baseline(1),0));
            XYp = R(-angle(obj.baseline(1),0)) * [X'; Y']; XYp(:,XYp(2,:)>centerp(2)) = [];
            Xp = XYp(1,:); Yp = XYp(2,:);
            
            %---------------------------------------------------------------
            %         Step 4 Separate left and right side and calculate
            %         average distance from center for both sides
            %---------------------------------------------------------------
            Xl = Xp(Xp < centerp(1)); Yl = Yp(Xp < centerp(1));
            Xr = Xp(Xp > centerp(1)); Yr = Yp(Xp > centerp(1));
            [Yl, order] = sort(abs(Yl - centerp(2))); Xl = centerp(1) - Xl(order);
            [Yr, order] = sort(abs(Yr - centerp(2))); Xr = Xr(order) - centerp(1);
            %Xc = mean([Xl; Xr],1); Yc = mean([Yl; Yr],1);
            
            %---------------------------------------------------------------
            %         Step 5 Numerically integrate cylinders of  distance 
            %                from center for both sides
            %---------------------------------------------------------------
            volume = trapz(Yl,pi*Xl.^2)/2 + trapz(Yr,pi*Xr.^2)/2;
        end
    end
    methods (Hidden = true)
        function quality = gaugeQuality(obj) % makes image corrections to improve accuracy and determines image quality
            %---------------------------------------------------------------
            %         Quality 1 Image passes enclosed edge criterion
            %---------------------------------------------------------------
            H = outlineDrop(obj); % conversion to grayscale, edge detection, fills in small gaps in edge, removes small false positives, adds bottom line
            I = imfill(logical(H),[1,1]); 	% fill image from top right
            
            if sum(I(:)) < 0.9*numel(I) % check if center region remains black else continue
                quality = 1; return
            end
            
            %---------------------------------------------------------------
            %         Quality 2 Droplet is darker than background
            %---------------------------------------------------------------
            %---------------------------------------------------------------
            %         Step 1 Smoothen image for image processing
            %---------------------------------------------------------------
            smoothIm = imopen(rgb2gray(obj.Image),strel('disk',15));
            adjIm = imsharpen(imadjust(smoothIm));
            
            %---------------------------------------------------------------
            %         Step 2 Evaluate image quality
            %---------------------------------------------------------------
            % find droplet edge and estimate approximate drop size by convex hull
            dropAprox = bwconvhull(edge(adjIm));
            % approximate background
            th = 8; % set threshold for background
            bg = grayconnected(adjIm,1,1,th) + grayconnected(adjIm,size(adjIm,1),1,th) + grayconnected(adjIm,1,size(adjIm,2),th) + grayconnected(adjIm,size(adjIm,1),size(adjIm,2),th);
            bg = bg > 0;
            % compare region brightness: good quality (drop darker bg, back illumination),
            % poor quality (drop brighter bg, front illumination)
            if mean(adjIm(dropAprox)) < mean(adjIm(bg))
                quality = 2;
            else
                quality = 3;
            end
        end
        function contactAngle = contactAngle(obj, ell, side) % calculates contact angles for a circle with baseline
            % calculates left and right contact angles from shape of
            % droplet
            
            %% ensure that a contact angle can be calculated and output is readable
            contactAngle = struct('CA',nan,'pin',[nan nan]);
            if ~isstruct(ell) || any( structfun(@isempty, ell) )
                return; 
            else
            end

            % angle returns the angle between two slopes m1, m2 in degrees
            angle = @(m1,m2)180-rad2deg(rem(pi+atan((m1-m2)/(1+m1*m2)),pi));
            
            %% equations to calculate sloped line intersections with a rotated ellipse
            % from
            % http://quickcalcbasic.com/ellipse%20line%20intersection.pdf
            % change from deg to rad?
            A = @(h,v,m,phi)v^2*(cosd(phi)^2 + 2*m*cosd(phi)*sind(phi) + m^2*sind(phi)^2) + ...
                            h^2*(m^2*cosd(phi)^2 - 2*m*cosd(phi)*sind(phi) + sind(phi)^2);
            B = @(h,v,m,b,phi)2*v^2*b*(cosd(phi)*sind(phi) + m*sind(phi)^2) + ...
                              2*h^2*b*(m*cosd(phi)^2 - cosd(phi)*sind(phi));
            C = @(h,v,b,phi)b^2*(v^2*sind(phi)^2 + h^2*cosd(phi)^2) - h^2*v^2;
            
            X = @(h,v,m,b,phi,sign) (-B(h,v,m,b,phi)  +sign*sqrt(B(h,v,m,b,phi)^2 - 4*A(h,v,m,phi)*C(h,v,b,phi))) / ...
                                        (2*A(h,v,m,phi));
            Y = @(x)obj.baseline(1)*x + obj.baseline(2);
            
            % from http://quickcalcbasic.com/ellipse%20line%20intersection.pdf
            m1 = @(h,v,e,f,x,y,phi)-((x-e)*cosd(phi) + (y-f)*sind(phi)) / ...
                                    ((y-f)*cosd(phi) - (x-e)*sind(phi)) * v^2/h^2;
            m2 = @(m1,phi)(m1*cosd(phi)+sind(phi))/(cosd(phi)-m1*sind(phi));
            m2p = @(phi)-cosd(phi)/sind(phi);
            
            
            %% contact angle
            % calculate intersection with baseline
            x = X(ell.a, ell.b, obj.baseline(1), ...
                  obj.baseline(2) + obj.baseline(1)*ell.X0_in - ell.Y0_in, ... % shited line intersection with Y-axis
                  -ell.phi, side) + ell.X0_in; % side -1 for left, +1 for right
            y = Y(x);
            contactAngle.pin = [x, y];
            
            % calculate tangent with baseline
            tangent = m2(m1(ell.a,ell.b,ell.X0_in,ell.Y0_in,x,y,-ell.phi),-ell.phi);
            tangent(isnan(tangent)) = m2p(ell.phi);
            tangent(~isfinite(tangent)) = sign(tangent)*1e99; % replacing infinity
            
            if isreal(tangent)
                if side < 0
                    contactAngle.CA = angle(tangent,obj.baseline(1));
                else
                    contactAngle.CA = angle(obj.baseline(1),tangent);
                end
            end
            
%             [x, y] =  linecirc(obj.baseline(1),obj.baseline(2),circle.center(1),circle.center(2), circle.radius); % intersections circle with baseline
%             for n = 1:length(x)
%                 tangentDirection = [x(n)-circle.center(1), y(n)-circle.center(2)]*[0 -1; 1 0]; % tangent direction
%                 tangent = tangentDirection(2)/tangentDirection(1);
%                 tangent(~isfinite(tangent)) = sign(tangent)*1e99; % replacing infinity
% 
%                 if x(n) < circle.center(1)
%                     CAs(1) = angle(tangent,obj.baseline(1));
%                 else
%                     CAs(2) = angle(obj.baseline(1),tangent);
%                 end
%                 contactAngle = CAs(side);
%             end
        end 
        function H = outlineDrop(obj)
            E = edge(rgb2gray(obj.Image));  % classic edge detection in gray-scale image
            F = bwmorph(E,'bridge',Inf);    % connect small gaps in edge detected 
            G = bwareaopen(F,25);           % remove small areas
            line = @(x)round(obj.baseline(2) + obj.baseline(1)*x);
            H = G; H(sub2ind(size(H),line(1:size(E,1)), 1:size(E,1))) = 1;      
            % H(size(E,1),:) = 1;      % add bottom white line to close edge
        end
        function [Pl, Pr] = pointsFurthestApart(obj)
            angle = @(m1,m2)180-rad2deg(rem(pi+atan((m1-m2)/(1+m1*m2)),pi));
            R = @(phi)[ cosd(phi) sind(phi); -sind(phi) cosd(phi) ]; % rotation matrix in 2D
            X = @(a,X0,range)  a * cos(range) + X0; Y = @(b,Y0,range) b * sin(range) + Y0; % functions plotting X,Y-values of ellipse
            rotEll = @(ell,range) R(ell.phi) * [X(ell.a,ell.X0,range); Y(ell.b,ell.Y0,range)];
            
            try E1 = R(angle(0,obj.baseline(1))) * rotEll(obj.Shape(1),(0 : 0.01 : 2*pi) + angle(0,obj.baseline(1)));
                E2 = R(angle(0,obj.baseline(1))) * rotEll(obj.Shape(3),( 0 : 0.01 :  2*pi ) + angle(0,obj.baseline(1)));
                
                % cut left and right edge of ellipse
                E1 = E1(:, all( E1 < (R(angle(0,obj.baseline(1)))*obj.CA.leftPin' + [10; 0]), 1));
                E2 = E2(:, arrayfun(@(row) all(E2(:,row) < (R(angle(0,obj.baseline(1)))*obj.CA.rightPin' - [10; 0]) == [0; 1]), 1:size(E2,2)) );
                
                % interpolate values to set y-values
                query = linspace(max([min(E1(2,:)) min(E2(2,:))]), min([max(E1(2,:)) max(E2(2,:))]), 200);
                E1i = [interp1(E1(2,:),E1(1,:),query); query];
                E2i = [interp1(E2(2,:),E2(1,:),query); query];
                    
                [~, ind] = max(abs(E2i(1,:) - E1i(1,:)));
                
                Pl = R(-angle(0,obj.baseline(1)))*E1i(:,ind);
                Pr = R(-angle(0,obj.baseline(1)))*E2i(:,ind);
            catch
                Pl = [NaN NaN]; Pr = [NaN NaN];
            end
        end
    end
    methods (Static, Hidden = true)
        function dist = maxDistanceBetween(circle1, circle2) % calculates outer distance of two circles
            try minX = circle1.center(1)-circle1.radius;
                maxX = circle2.center(1)+circle2.radius;
                dist = abs(maxX-minX)/2;
            catch
                dist = NaN;
            end
        end
        function ellipse = createEllipse(radii,center,phi)
            try ellipse.a     = radii(1);   % sub axis (radius) of the X axis of the non-tilt ellipse
                ellipse.b     = radii(end); % sub axis (radius) of the Y axis of the non-tilt ellipse
                ellipse.X0    = center(1);  % center at the X axis of the non-tilt ellipse
                ellipse.Y0    = center(2);  % center at the Y axis of the non-tilt ellipse
                ellipse.phi   = rad2deg(phi);        % orientation in radians of the ellipse (tilt)
                
                R             = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
                P_in          = R * [center(1);center(2)];
                ellipse.X0_in = P_in(1);    % center at the X axis of the tilted ellipse
                ellipse.Y0_in = P_in(2);    % center at the Y axis of the tilted ellipse
            catch
                ellipse.a     = nan;
                ellipse.b     = nan;
                ellipse.X0    = nan;
                ellipse.Y0    = nan;
                ellipse.phi   = nan;
                ellipse.X0_in = nan;
                ellipse.Y0_in = nan;
            end
        end
        function [centers, radii, metric] = detectCircles(Im,radius,varargin)
            % finds circles in Im using Matlab Hough transform algorithm

            % validate input
            validateattributes(Im, {'numeric', 'uint8'},{'3d'})
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
    end
end