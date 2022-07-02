function project(~)
clc

I1 = imread("BIM472_Image12.jpg");

xdim = size(I1, 2);
ydim = size(I1, 1);

R1 = zeros(ydim, xdim, 3, 'uint8');

% Image with only red channel
T1 = zeros(ydim, xdim, 1);

% Blank image where polygon is drawn
matchData = zeros(2, 4); 

% Store information on overlap
hsvI = rgb2hsv(I1);

N_sides = 4; 
% Number of sides in polygon
t = (1 /(N_sides * 2):1 / N_sides:1)' * 2 * pi; 
% t vector for polygon
iteration = 1; 
% number of times images has been checked for polygon fit

HueLower = 20 / 255;
HueUpper = 150 / 255;
SaturationMin = 50 / 255;
ValueMin =  50 / 255;


for x = 1:1:xdim
    for y = 1:1:ydim
        if (hsvI(y, x, 1) > HueUpper || hsvI(y, x, 1) < HueLower) && (hsvI(y, x, 2) > SaturationMin) && (hsvI(y, x, 3) > ValueMin)
            R1(y, x, 1) = I1(y, x, 1);
            R1(y, x, 2) = I1(y, x, 2);
            R1(y, x, 3) = I1(y, x, 3);
        end
    end
end

E1 = edge(rgb2gray(R1), 'sobel', 0.1);

for x = 2:1:xdim
    for y = 2:1:ydim
        if E1(y, x) == 1
           T1(y + 1, x - 1) = 1;
           T1(y + 0, x - 1) = 1;
           T1(y - 1, x - 1) = 1;
           T1(y + 1, x + 0) = 1;
           T1(y - 0, x + 0) = 1;
           T1(y + 1, x + 1) = 1;
           T1(y + 0, x + 1) = 1;
           T1(y - 1, x + 1) = 1;
        end
    end
end

s = strel('square', 5);
C1 = imclose(T1, s);
BW2 = imfill(C1,'holes');

% Get largest objects
F1 = bwareafilt(logical(BW2), 3); 
% use 8 largest objects as safety margin
N1 = bwareafilt(logical(T1), 3); 
% Image with blank pixels removed

% Get boundaries
[~, ~, ~, A] = bwboundaries(E1, 'noholes'); 
% B = Boundaries, N = number of objects, 4 = connectivity

% Only grab slices of images that contain pixels to improve efficiency
columns = any(N1);
Loc = 1;

tic
for x = 1:1:xdim
        if columns(x) == 1
            Slice1(:, Loc) = N1(:, x);
            Slice1Image(:, Loc, 1) = I1(:, x, 1);
            Slice1Image(:, Loc, 2) = I1(:, x, 2);
            Slice1Image(:, Loc, 3) = I1(:, x, 3);
            Loc = Loc + 1;
        end
end
toc

rows = any(Slice1, 2);
Loc = 1;

tic
for y = 1:1:ydim
        if rows(y) == 1  
            Slice2(Loc, :) = Slice1(y, :);
            Slice2Image(Loc, :, 1) = Slice1Image(y, :, 1);
            Slice2Image(Loc, :, 2) = Slice1Image(y, :, 2);
            Slice2Image(Loc, :, 3) = Slice1Image(y, :, 3);
            Loc = Loc + 1;
        end
end
toc

% Reduce resolution to further improve efficiency
multiplier = 128 / size(Slice2, 1);
if multiplier < 1
    Slice2 = imresize(Slice2, multiplier);
end

xdim = size(Slice2, 2);
ydim = size(Slice2, 1);

% Cycle through polygon scale factors and position to determine if there's a good match
for scale = (xdim / 5):3:(xdim)
    for xOffset = abs(min(scale*sin(t))):3:xdim - abs(min(scale*sin(t)))
        for yOffset = abs(min(scale*cos(t)))+ 1:3:ydim - max(scale*cos(t))         
            % Get vertices
            polyCoord = [round(scale*sin(t) + xOffset) round(scale*cos(t) + yOffset)];
            polyCoord(N_sides + 1, 1) = polyCoord(1, 1);
            polyCoord(N_sides + 1, 2) = polyCoord(1, 2);
            
            % Create polygon using vertices
            j = 1;
            k = 1;
            
            for i = 1:1:(N_sides * 2) 
                if mod(i, 2) == 1 
                    poly(i) = polyCoord(j, 1);
                    j = j + 1;
                else
                    poly(i) = polyCoord(k, 2);
                    k = k + 1;
                end
            end
            
            % Draw Polygon on fresh image 
            S1 = zeros(ydim, xdim, 3, 'uint8');
            S1 = insertShape(S1, 'Polygon', poly, 'Color', 'white');

            % Convert image to logical type
            S2 =  imbinarize(rgb2gray(S1));

            % Get overlap
            Slice2 = imcrop(Slice2,[0 0 xdim ydim]);
            M1 = Slice2 .* S2;

            % Get number of overlapped pixels in images
            commonPixels = sum(M1(:) == 1);
            
            % Save data
            matchData(iteration, :) = [commonPixels scale xOffset yOffset];
            iteration = iteration + 1;     
        end
    end
end

% Calculate best match
[~, I] = max(matchData);
scale = matchData(I(1), 2);
xOffset = matchData(I(1), 3);
yOffset = matchData(I(1), 4);

% Create polygon for Best Fit Shape overlapped on Sliced/Processed Image
% Get vertices
polyCoord = [round(scale*sin(t) + xOffset) round(scale*cos(t) + yOffset)];
polyCoord(N_sides + 1, 1) = polyCoord(1, 1);
polyCoord(N_sides + 1, 2) = polyCoord(1, 2);

% Create polygon using vertices
j = 1;
k = 1;

for i = 1:1:(N_sides * 2) 
    if mod(i, 2) == 1 % odd
        poly(i) = polyCoord(j, 1);
        j = j + 1;
    else
        poly(i) = polyCoord(k, 2);
        k = k + 1;
    end
end

% Create another polygon for Best Fit Shape overlapped on Sliced Raw Image
% Only needed if multiplier < 1
if multiplier < 1
    polyCoord2 = [round((scale*sin(t) + xOffset)/multiplier) round((scale*cos(t) + yOffset)/multiplier)];
    polyCoord2(N_sides + 1, 1) = polyCoord2(1, 1);
    polyCoord2(N_sides + 1, 2) = polyCoord2(1, 2);
else
    polyCoord2 = polyCoord;
end

% Create polygon using vertices
j = 1;
k = 1;

for i = 1:1:(N_sides * 2) 
    if mod(i, 2) == 1 % odd
        poly(i) = polyCoord2(j, 1);
        j = j + 1;
    else
        poly(i) = polyCoord2(k, 2);
        k = k + 1;
    end
end

sum(A(:) == 1);

figure('Name',sprintf('BIM472_Image%02d.jpg', 1),'NumberTitle','off');

subplot(2, 3, 1), imshow(I1);
title('Original');

subplot(2, 3, 2), imshow(R1);
title('Extract Red Pixels');


subplot(2, 3, 3), imshow(T1);
title('Fill in Edges');

subplot(2, 3, 4), imshow(F1);
title('Largest Object');

subplot(2, 3, 5)
hold on;
imshow(Slice2Image);
plot(polyCoord2(:,1), polyCoord2(:,2), 'g')
hold off;
title('Best Fit Shape');
