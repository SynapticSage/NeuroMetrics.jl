%% Goal Coding in CA1
%% Analyze
% Purpose of this script is to run a filter-framework driven exploration of 
% all CA1 cells, per animal, and examine tuning curves to goal distance, goal 
% linear distance, and goal angle.
% 
% Days 1 - 5 , Normal W - Track , Run Sessions
% 
% Cells in $CA1_{dorsal}$ and $CA1_{intermediate},$ whose mean firing rate 
% $f$ is $f >  \frac{1}{2} hz$

for animal = {'HPa','HPb','HPc'}

a = animal{1};
folder=['~/Data/Local/Current/GoalCoding/' 'PFC_' a];
mkdir(folder);
pushd(folder);

days = '1:5';
epochs='isequal($environment,''wtr1'') && isequal($type,''run'')';
% cells = ['(isequal($area,''CA1'') || ' ...
%          'isequal($area,''iCA1'') ) && '...
%          '($meanrate > 0.5)'];
cells = ['isequal($area,''PFC'') &&'...
         '($meanrate > 0.5)'];

%% 
% Now we instantiate our filter to run these characteristics,

modf = createfilter('animal',animal,'days',days,'epochs',epochs,'cells',cells,'iterator','singlecellanal');
modf=setfilterfunction(modf,'DFA_spikingtrigger_v2',{'pos','linpos','spikes','cellinfo'},30e-3,'ploton',true,'mode','separate');
modf= runfilter(modf);
%% Pre-plot
% (1) Set home directory for results
%% 
% (2) Save outputs

extract_modf(modf,'mode','store_nosarray','append',{'animal',a});
popd;
end
