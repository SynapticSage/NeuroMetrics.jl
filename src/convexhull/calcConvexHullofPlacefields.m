function [ out ] = calcConvexHullofPlacefields( mapfield, varargin )
%PLACEFIELD Allows one to ID the times an animal
%is is some subset of placefields. It can be used either as a time filter
%for all cells or as a stand-alone function to store per all cells times
%when the animal is inside the convex hull defining a field
%
%   Per listed cell one or more convex hulls of the place field(s) are
%   found. It does this under control of a threshold parameter: either raw
%   value or normalized value. Then. this algorithm finds the set of convex
%   hulls for that threshold. If it's more than one, a condition lets you
%   filter out the largest convex hull.
%
%   --- Params ---
%   
% TODO Sarel et al 2017 -- add expansion of hull by 2 bins to ensure no
% in-field spikes

threshold = 1; %hz/occ
thresholdtype = 'halfpeak';
largestneighborhoodnum = inf; % if a number less than inf, algo finds the largest n connected regions, before finding the convex hull
optlistassign(who,varargin);

% Take the mapfields and acquire occ norm firing rate
% (1) Acquire the spike rate
sr = mapfield.smoothedspikerate;
% (2) Occupancy
occ = mapfield.smoothedoccupancy;
% (3) Occ norm fields
occsr = sr ./ occ;
% (4) Record where values are negative or infinite and remove them.
impossiblevals = occsr==inf|occsr==-inf|occsr<0;
occsr(impossiblevals) = 0;

if strcmp(thresholdtype,'halfpeak') % same used by Sarel et al. 2017
    threshold = max(max(mapfield.spikerate))/2;
%       threshold = max(max(occsr))/2;
end

% Now, figure all above threshold values
abovethreshold = occsr >= threshold;
abovethreshold = mapfield.spikerate >= threshold;

% Filter out disconnected regions
if largestneighborhoodnum<inf
    % Get list of contiguous regions
    areas = bwconncomp(abovethreshold);
    % Find the N largest
    areasizes = cellfun(@numel,areas.PixelIdxList);
    % Which largest?
    areasizes = [1:numel(areasizes); areasizes];
    areasizes=sortrows(areasizes',2);
    largestregions = areasizes(1:largestneighborhoodnum,1);
    % Now zero out anything not apart of these large regions
    abovethreshold(cat(1,areas.PixelIdxList{largestregions}))=0;
end
[aty, atx]= find(abovethreshold);


out.occspikerate = occsr;
out.spikerate = mapfield.smoothedspikerate;
out.abovethresh = abovethreshold;
out.xticks=mapfield.xticks;
out.yticks=mapfield.yticks;

try
% Last, find the convex hull that encapsulates the set of all above
% threshold points
hull = convhull(atx,aty);
out.indconvhull = [atx(hull), aty(hull)];

% Now just to be nice, lets include a polygon of indices AND of the xy
% animal positions, so that users can easily find points in the place
% field using polygonin
out.posconvhull = [mapfield.xticks(atx(hull)); mapfield.yticks(aty(hull))]';
catch
  fprintf('cannot find hull\n');
  out.indconvhull=[];
  out.posconvhull=[];
end

out.settings.largestneighborhoodnum = largestneighborhoodnum;
out.settings.threshold = threshold;
out.settings.thresholdtype = thresholdtype;
out.settings.dateproduced = date;
