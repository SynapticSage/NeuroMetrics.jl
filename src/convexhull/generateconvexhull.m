% generateconvhull.m
% Generates convex hull for a given animal
function generateconvexhull(animalprefix,days,varargin)

  ainfo = animaldef(animalprefix);
  mapfields = loaddatastruct(ainfo{2:3},'mapfields');
  searchresults = cellfetch(mapfields,'xticks');
  if nargin >= 2
    badrows = ~ismember(searchresults.index(:,1),days);
    searchresults.index(badrows,:)=[]; % delete unrequested indices
  end  
  out=calcConvexHullofPlacefields(ainfo{2:3},searchresults.index,'mapfields',mapfields,'largestneighborhoodnum',1,varargin{:});
  
  savedatastruct(out,ainfo{2:3},'convexhull');
  
end