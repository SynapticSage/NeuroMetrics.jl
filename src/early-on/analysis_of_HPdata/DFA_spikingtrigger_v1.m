function out = DFA_spikingtrigger_v1(index,exclude,pos,linpos,spikes,cellinfo,varargin)

% TODO
% ------------------------
% Expand to track per well
% pool gangle versus rate
% pool gdist versus rate
% pool glindist versus rate
%
% pool within each trial, then across trials
% normalize by the number of times an animal is located at that state
% 
% May have added these, but now need to test that the added
% functionality works!

%% Parse parameters
ip = inputParser;
ip.addParameter('window',[0 10]); % controls windows to pull!
ip.addParameter('goalanal',true); 
ip.addParameter('mode','together'); % Goals 'together' or goals 'separate' , in separate, each well is treated as a separate goal when pulling triggers
ip.addParameter('ploton',false); 
ip.addParameter('nbin',100); % How many bins to use for the plot mode
ip.parse(varargin{2:end});

goalanal = ip.Results.goalanal;
window = ip.Results.window;
mode = ip.Results.mode;
ploton = ip.Results.ploton;
nbin = ip.Results.nbin;

cellinfo = cellinfo{index(1,1)}{index(1,2)}{index(1,3)}{index(1,4)};

%% Pos

% Extract pos variables
pos = pos{index(1,1)}{index(1,2)};
linpos = linpos{index(1,1)}{index(1,2)};
post = pos.data(:,1);
posd = pos.data(:,colof('dir',pos));
try
    posx = pos.data(:,colof('x-sm',pos));
    posy = pos.data(:,colof('y-sm',pos));
catch
    posx = pos.data(:,colof('x',pos));
    posy = pos.data(:,colof('y',pos));
end

% Calculate relevent goal variables for all times
angle_goal  = angle2goal(posx,posy,posd,linpos);
dist_goal   = dist2goal(posx,posy,linpos);
ldist_goal  = linearized2goal(linpos);

%% Firing

% Now determine measures for all times, given Ws found above
[F,time_F,T] = ry_calcrate(index, [], spikes, varargin{1}); % MAKE SURE THIS IS PROPERLY EXCLUDING!
S = F.*T; % Spike count
% F = (F-nanmean(F,2))./nanstd(F,1,2); % Normalize firing rate

%% Get Periods
% Use time filtered exclude times to separate
goodtimes_F = ~isExcluded(time_F,exclude);
goodtimes_pos = ~isExcluded(post,exclude);

% From here, iterate each inclusion period and then use interp to ensure
% each matches in its total duration
[inclusion_period_well,trigger_time] = getWellPeriod(post,linpos,window,mode);

%% Iterate and collect windows

% Windows times to interp
win_pos = -window(1):median( diff(post) ):window(2);
win_F = -window(1):median( diff(time_F) ):window(2);

windows = struct('F',[],'x',[],'y',[],'angle_goal',[],'dist_goal',[],'ldist_goal',[]);
windows = repmat(windows,numel(trigger_time),1);
for jj = 1:numel(trigger_time) % For each well/goal type
    for ii = 1:size(inclusion_period_well{jj},1) % For each time the animal visits that well, capture a window
        
        period = inclusion_period_well{jj}(ii,:);
        trigger = trigger_time{jj}(ii);

        % Obtain time slices
        idxf = goodtimes_F' & ...
            ( time_F >= period(1) & time_F <= period(2) );
        idxp = goodtimes_pos & ...
            ( post >= period(1) & post <= period(2) );

        % Interp a window of firing
        windows(jj).F = [windows(jj).F ; ... RATE
                interp1(time_F(idxf)-trigger,F(idxf),win_F) ];

        % Interp windows of x,y,and goal variables
        windows(jj).x = [windows(jj).x ; ...
                interp1(post(idxp)-trigger,posx(idxp),win_F) ];
        windows(jj).y = [windows(jj).y ; ...
                interp1(post(idxp)-trigger,posy(idxp),win_F) ];
        windows(jj).angle_goal = [windows(jj).angle_goal ; ...
                interp1(post(idxp)-trigger,angle_goal(idxp),win_F) ];
        windows(jj).dist_goal = [windows(jj).dist_goal ; ...
            interp1(post(idxp)-trigger,dist_goal(idxp),win_F) ];
        windows(jj).ldist_goal = [windows(jj).ldist_goal ; ...
            interp1(post(idxp)-trigger,ldist_goal(idxp),win_F) ];
    end
end

%% Output
% Store firing rate and independent principal information
out.F = F;
out.t = time_F;

% Store X,Y,Vel components
out.windows=windows;

out.day = index(1);
out.epoch = index(2);

if ploton || analyze
    
    ax = analyze_window(windows,nbin);
    
end

%% --- Helper Functions ---
% ------------------------------
% -----------------------------

    function ax = analyze_window(windows,nbin)
        % Carries out analysis on each goal-triggered window
        
        %tuning_fit = @(phi,w,amp,baseline) amp*cos(w*theta + ...
        %   phi)-baseline;
        
        figure(1);clf;df =get(0,'defaultAxesFontSize') ;set(0,'defaultAxesFontSize',16);
        
        nWin=numel(windows);
        nAnal = 4;
        
        ax = tight_subplot(nAnal,nWin,[0.15,0.07],0.1,0.1);
        ax = reshape(ax, nWin , nAnal);
        
        % Time Cell?
        for i = 1:nWin
            
                % Time cell?
                a = ax(i,1);
                plot(a, win_F, sum(windows(i).F,1)/size(windows(i).F,1) ,'bo-' ); 
                title(a,'Time');
                %fminsearch
                
                % Angle cell?
                a = ax(i,2);
                [y,x] = group(windows(i).angle_goal,windows(i).F,nbin); 
                plot(a,x,y,'ro-');title(a,'Angle');
                set(a, 'xtick',[0 pi/2 pi 3*pi/2 2*pi]);
                set(a, 'xticklabel',{'0' '\pi/2' '\pi' '3\pi/2' '2\pi' });
                
                % Distance cell?
                a = ax(i,3);
                [y,x] = group(windows(i).dist_goal,windows(i).F,nbin); 
                plot(a,x,y,'go-');title(a,'Distance');

                % Linear distance cell?
                a = ax(i,4);
                [y,x] = group(windows(i).ldist_goal,windows(i).F,nbin);
                plot(a,x,y,'ko-'); title(a,sprintf('Trajectory Location'),'fontsize',14);
        end
        
        
        linkaxes( ax(:,1) ); linkaxes( ax(:,2) ); linkaxes( ax(:,3) ); linkaxes( ax(:,4) );
        set(ax(1,1),'xlim',[-inf inf]);set(ax(2,:),'xlim',[-inf inf]);set(ax(3,1),'xlim',[-inf inf]);set(ax(3,1),'xlim',[-inf inf]);
        arrayfun( @(x) ylabel(ax(1,x),'$\bar{r}$','interpreter','latex'), 1:size(ax,2));
        xlabel(ax(1,end),'well 1');xlabel(ax(2,end),'well 2');xlabel(ax(3,end),'well 3');
        
        st = sprintf('%s Day %d Ep %d Tet %d Cell %d',cellinfo.area,index(:));
        
        axis( suptitle(st), 'off');
        mkpushd([pwd filesep 'cell_summary']);
        saveThis(gcf,pwd,st,{'png','fig'});
        popd;
        set(0,'defaultAxesFontSize',df);
                
    end

    function [yint,center,ynint] = group(x,ynint,nbin)
        % group "groups" a set of properties by another set of properties.
        
        x = x(:); ynint=ynint(:);
        
        [N,edges,bin] = histcounts(x,nbin);
        for u = unique(bin)'
            if u == 0; continue; end;
            Y{u} = sum(ynint(bin==u),1)/sum(bin==u);
            % TODO would be easy to add a bootstrap here.
        end
        
        center = mean([edges(1:end-1);edges(2:end)],1);
        ynint = cat(1,Y{:});
        yint=interp1(center(~isnan(ynint)),ynint(~isnan(ynint)),center);
        
        % Bootstrap
        
    end

    function [trigger_wins,triggers] = getWellPeriod(t,linpos,window,mode)
        % Gets set of triggered window for each goal
        
        if nargin < 4
            mode='together';
        end
        
        % Use start/stop of well trajectories
        wellex = linpos.statematrix.wellExitEnter(:,2);
        wellex_change = [0; diff(wellex)];
        
        % Acquire trigger times either per well or for all wells
        % collectively
        triggers = {};
        switch mode
            case {'all','together'}
            triggers{1} = wellex_change ~= 0 ;
            case 'separate'
            unique_well = sort(unique(wellex));
            for u = unique_well'
                triggers{end+1} = u == wellex & wellex_change ~= 0;
            end
            otherwise, error('Bad input');
        end
        
        for i = 1:numel(triggers)
            trigger_wins{i} = [t(triggers{i})-window(1), t(triggers{i})+window(2)];
            triggers{i} = t(triggers{i});
        end
        
    end

    % --------------------------------------------------------------------
    % Goal Distance/Angle/Linear Helper function
    % --------------------------------------------------------------------
    function dist = dist2goal(x,y,linpos)
        % Calculates the non-linearized distance to all three wells
        
        well=linpos.wellSegmentInfo.wellCoord;
        
        x_well = well(:,1);
        y_well = well(:,2);
        
        dist=sqrt( bsxfun(@minus,x,x_well').^2 + bsxfun(@minus,y,y_well').^2 );
        
    end

    function linearized = linearized2goal(linpos)
        % Calculates the linearized distance to all three wells
        linearized=linpos.statematrix.linearDistanceToWells;
    end
    
    function angle = angle2goal(x,y,dir,linpos)
        % Calculates the angle between the animal and all of the wells
        
        % Get well coordinates
        well=linpos.wellSegmentInfo.wellCoord;
        % Acquire HEADING direction
        x_heading = [0; diff(x)];
        y_heading = [0; diff(y)];
        % Heading dir
        dir_heading = atan2(y_heading,x_heading);
        
        % Get x and y of each well
        x_well = well(:,1);
        y_well = well(:,2);
 
        % Calculate angle
        angle_well = atan2( ...
            bsxfun(@minus,y_well',y),...
            bsxfun(@minus,x_well',x)...
            );
        
%         [dir_u, dir_v]= unitvector(dir);
%         [angle_well_u, angle_well_v]= unitvector(angle_well);
%         
%         vec{1} = [dir_u,dir_v];
%         vec{2} = [angle_well_u(:,1), angle_well_v(:,1)];
%         
%         plot2dVector(x,y,vec,8)
        
        % Calculate angle to well
        angle = quad2normal(bsxfun(@minus,dir_heading,angle_well));
        
    end

    function quad = quad2normal(quad)
        quad(quad<0) = quad(quad<0)+2*pi;
    end

    function [u,v] = unitvector(dir)
        
        dir = ( (dir-min(dir(:)))/abs(range(dir(:))) * 2*pi );
        u=cos(dir);
        v=sin(dir);
        
    end

    % --------------------------------------------------------------------
    % Time period helper function
    % --------------------------------------------------------------------
    function includePeriods = getInclusionRange(time,included)
       % includePeriods = getInlcudePeriods(time, included)
        % Calculates the start and end times for all exclusion periods, given a
        % time vector and an include vector of 1s and 0s of the same length.
        
        assert(isrow(included)||iscolumn(included),'Your control statement outputs more than a single row or column.');
        
        % Ensure time is a column
        if isrow(time), time=time(:); end

        % Check for equal length of included binary vector and time
        if (length(time) ~= length(included))
            error('The TIME and INCLUDED vectors must me the same length');
        end

        % Discover start and stop times
        starttimes = find((diff(included) == 1))+1;
        starttimes = starttimes(:);
        endtimes = find((diff(included) == -1));
        endtimes = endtimes(:);
        if (included(1) == 1)
            starttimes = [1; starttimes];
        end
        if (included(end) == 1)
            endtimes = [endtimes; length(included)];
        end
        
        if numel(starttimes) > numel(endtimes)
            if starttimes(1) < endtimes(1)
                starttimes(1) = [];
            else
                starttimes(end) = [];
            end
        elseif numel(endtimes) > numel(starttimes)
            if endtimes(1) < starttimes(1)
                endtimes(1) = [];
            else
                endtimes(end) = [];
            end
        end
    
        % Create inclusion vector
        includePeriods = [time(starttimes) time(endtimes)]; 
        assert( sum(starttimes-endtimes) <= 0 );
        
    end

end
