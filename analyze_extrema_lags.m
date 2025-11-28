function analyze_extrema_lags()
%% ================= USER INPUTS ==========================================
xlsxFile   = 'Results_30min_60min.xlsx';   % file name
sheetNames = {'30min','60min'};           % first 30, then 60 for subplots

% Which discharge series to use as "Qmain" (e.g., 'Obs_Discharge', 'Benchmark', 'Model')
dischargeField = 'ObsQ';

% Warmup period (days): ignore data earlier than this
warmupDays = 14;

% Peak detection sensitivity (adjust if you get too many / too few extrema)
prominenceFrac = 0.05;   % as fraction of data range

%% === GLOBAL PLOT STYLE (Montserrat + LaTeX, ticks out, thick borders) ===
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% Requires Montserrat installed on your system
set(groot,'defaultAxesFontName','Montserrat');
set(groot,'defaultTextFontName','Montserrat');

%% Colors: blue (Qmain), cyan (Qmodel), red (Ksat), dark gray (Sy)
colQmain  = [0 0 1];        % blue
colQmodel = [0 1 1];        % cyan
colKsat   = [1 0 0];        % red
colSy     = [0.3 0.3 0.3];  % dark gray

%% ========================================================================
figure('Name','Normalized Q, Q_{model}, K_{sat}, S_y extrema (30 vs 60 min)', ...
       'Color','w');

for s = 1:numel(sheetNames)
    sheet = sheetNames{s};
    fprintf('\n===== Sheet: %s =====\n',sheet);

    % Read sheet as table
    T = readtable(xlsxFile,'Sheet',sheet);

    % Variable names (readtable replaces spaces with underscores)
    fprintf('Variables in this sheet:\n');
    disp(T.Properties.VariableNames);

    % Extract series (time in days)
    tAll      = T.Time;
    qMainAll  = T.(dischargeField);  % main discharge series (e.g., obs/benchmark)
    qModelAll = T.Qmodel;            % modeled discharge column
    ksatAll   = T.Ksat;
    syAll     = T.Sy;

    % ---- Apply warmup: use only data AFTER warmupDays ----
    mask   = tAll >= warmupDays;
    t      = tAll(mask);
    qMain  = qMainAll(mask);
    qModel = qModelAll(mask);
    ksat   = ksatAll(mask);
    sy     = syAll(mask);

    if numel(t) < 5
        warning('After warmup there are very few points in %s; check time units.',sheet);
    end

    % ---- Compute KGE between Qmain (obs) and Qmodel (sim) ----
    KGE = calcKGE(qMain, qModel);
    fprintf('KGE (Qmain vs Qmodel) for %s: %.3f\n', sheet, KGE);

    % Choose prominence based on each series range
    promQmain  = (max(qMain)  - min(qMain))  * prominenceFrac;
    promQmodel = (max(qModel) - min(qModel)) * prominenceFrac;
    promKsat   = (max(ksat)   - min(ksat))   * prominenceFrac;
    promSy     = (max(sy)     - min(sy))     * prominenceFrac;

    % ---- Find peaks & valleys (raw values) ----
    extQmain  = findExtrema(t, qMain,  promQmain);
    extQmodel = findExtrema(t, qModel, promQmodel); %#ok<NASGU> % kept for possible later use
    extKsat   = findExtrema(t, ksat,   promKsat);
    extSy     = findExtrema(t, sy,     promSy);

    % ---- Compute lags (in days + hours) ----
    lagKsat_Qmain  = computeNearestLag(extQmain,  extKsat, 'Ksat');
    lagSy_Qmain    = computeNearestLag(extQmain,  extSy,   'Sy');
    lagKsat_Qmodel = computeNearestLag(findExtrema(t,qModel,promQmodel), extKsat, 'Ksat');
    lagSy_Qmodel   = computeNearestLag(findExtrema(t,qModel,promQmodel), extSy,   'Sy');

    % Basic pattern info (days)
    fprintf('Mean lag Qmain vs Ksat (hours):  %.3f\n', mean(lagKsat_Qmain.Lag_days*24));
    fprintf('Mean lag Qmain vs Sy   (hours):  %.3f\n', mean(lagSy_Qmain.Lag_days*24));
    fprintf('Mean lag Qmodel vs Ksat (hours): %.3f\n', mean(lagKsat_Qmodel.Lag_days*24));
    fprintf('Mean lag Qmodel vs Sy   (hours): %.3f\n', mean(lagSy_Qmodel.Lag_days*24));

    % Expose results in base workspace
    suffix = strrep(sheet,'-','_'); % safe for variable names
    assignin('base',['lagKsat_Qmain_'  suffix],lagKsat_Qmain);
    assignin('base',['lagSy_Qmain_'    suffix],lagSy_Qmain);
    assignin('base',['lagKsat_Qmodel_' suffix],lagKsat_Qmodel);
    assignin('base',['lagSy_Qmodel_'   suffix],lagSy_Qmodel);

    assignin('base',['extQmain_'  suffix],extQmain);
    assignin('base',['extKsat_'   suffix],extKsat);
    assignin('base',['extSy_'     suffix],extSy);

    %% ---- Normalization to [0,1] (per variable) ----
    [normQmain,  qMainMin,  qMainRange]   = normalize01(qMain);
    [normQmodel, qModelMin, qModelRange]  = normalize01(qModel);
    [normKsat,   ksatMin,   ksatRange]    = normalize01(ksat);
    [normSy,     syMin,     syRange]      = normalize01(sy);

    % Normalized values for extrema (for plotting)
    normQmain_ext  = (extQmain.Value  - qMainMin)  / qMainRange;
    normKsat_ext   = (extKsat.Value   - ksatMin)   / ksatRange;
    normSy_ext     = (extSy.Value     - syMin)     / syRange;

    %% ---- Plot in subplot (1: 30min, 2: 60min) ----
    subplot(2,1,s);

    % Lines
    plot(t, normQmain,  '-','Color',colQmain,  'LineWidth',1.8); hold on
    plot(t, normQmodel, '--','Color',colQmodel,'LineWidth',1.8);
    plot(t, normKsat,   '-.','Color',colKsat,  'LineWidth',1.8);
    plot(t, normSy,     ':','Color',colSy,     'LineWidth',1.8);

    % Peaks and valleys: Qmain
    plot(extQmain.Time(extQmain.Type=="peak"), ...
         normQmain_ext(extQmain.Type=="peak"), ...
         'o','MarkerFaceColor',colQmain,'MarkerEdgeColor','k','LineWidth',1);
    plot(extQmain.Time(extQmain.Type=="valley"), ...
         normQmain_ext(extQmain.Type=="valley"), ...
         'v','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colQmain,'LineWidth',1.2);

    % Peaks and valleys: Ksat
    plot(extKsat.Time(extKsat.Type=="peak"), ...
         normKsat_ext(extKsat.Type=="peak"), ...
         's','MarkerFaceColor',colKsat,'MarkerEdgeColor','k','LineWidth',1);
    plot(extKsat.Time(extKsat.Type=="valley"), ...
         normKsat_ext(extKsat.Type=="valley"), ...
         'd','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colKsat,'LineWidth',1.2);

    % Peaks and valleys: Sy
    plot(extSy.Time(extSy.Type=="peak"), ...
         normSy_ext(extSy.Type=="peak"), ...
         '^','MarkerFaceColor',colSy,'MarkerEdgeColor','k','LineWidth',1);
    plot(extSy.Time(extSy.Type=="valley"), ...
         normSy_ext(extSy.Type=="valley"), ...
         'v','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colSy,'LineWidth',1.2);

    ylabel('Normalized value (0-1)')
    if s == numel(sheetNames)
        xlabel('Time (days)','Interpreter','latex');
    else
        xlabel('','Interpreter','latex');
    end

    % Title with KGE
    titleStr = sprintf(['Normalized $Q$, $Q_{\\mathrm{model}}$, ', ...
                        '$K_{\\mathrm{sat}}$, $S_y$ extrema - %s ', ...
                        '(KGE = %.2f)'], sheet, KGE);
    title(titleStr,'Interpreter','latex');

    % Axes styling with daily grid only
    ax   = gca;
    xMin = floor(min(t));
    xMax = ceil(max(t));
    ax.YLim = [0 1];           % fixed for normalized
    styleAxesDaily(ax,xMin,xMax);

    if s == 1
        legend({'$Q$ (norm)', '$Q_{\mathrm{model}}$ (norm)', ...
                '$K_{\mathrm{sat}}$ (norm)', '$S_y$ (norm)', ...
                '$Q$ peaks', '$Q$ valleys', ...
                '$K_{\mathrm{sat}}$ peaks', '$K_{\mathrm{sat}}$ valleys', ...
                '$S_y$ peaks', '$S_y$ valleys'}, ...
                'Location','best','Interpreter','latex');
    end
end

end

%% ================= HELPER FUNCTIONS =====================================

function ext = findExtrema(t,y,prom)
% Find peaks and valleys of a time series y(t)
% t  : time vector (days)
% y  : series
% prom : MinPeakProminence

    % Peaks
    [pk, locPk] = findpeaks(y,t,'MinPeakProminence',prom);

    % Valleys -> peaks of -y
    [valNeg, locVal] = findpeaks(-y,t,'MinPeakProminence',prom);
    val = -valNeg;

    % Build table
    Time  = [locPk; locVal];      % days
    Value = [pk;    val];
    Type  = [repmat("peak",  numel(pk),1);
             repmat("valley",numel(val),1)];

    % Sort by time
    [Time,idx] = sort(Time);
    Value = Value(idx);
    Type  = Type(idx);

    ext = table(Time,Value,Type, ...
        'VariableNames',{'Time','Value','Type'});
end

function lagTable = computeNearestLag(extQ,extPar,parName)
% For each discharge extremum, find nearest extremum in parameter (Ksat/Sy)
% and compute time lag (days + hours).

    nQ = height(extQ);
    nearestTime   = zeros(nQ,1);
    nearestType   = strings(nQ,1);
    lagDays       = zeros(nQ,1);

    for i = 1:nQ
        t0 = extQ.Time(i);
        [~, idx] = min(abs(extPar.Time - t0));
        lagDays(i)      = extPar.Time(idx) - t0;  % +: parameter extremum AFTER Q
        nearestTime(i)  = extPar.Time(idx);
        nearestType(i)  = extPar.Type(idx);
    end

    lagHours = lagDays*24;

    lagTable = table( ...
        extQ.Time, extQ.Type, ...               % discharge extremum time & type
        nearestTime, nearestType, ...           % parameter extremum
        lagDays, lagHours, ...
        'VariableNames', { ...
            'Q_Time_days','Q_Type', ...
            [parName '_Time_days'],[parName '_Type'], ...
            'Lag_days','Lag_hours'});
end

function [yNorm, yMin, yRange] = normalize01(y)
% Normalize vector y to [0,1] using min–max.
    yMin   = min(y);
    yMax   = max(y);
    yRange = yMax - yMin;
    if yRange == 0
        yRange = 1;           % avoid division by zero; all values become 0
    end
    yNorm = (y - yMin) / yRange;
end

function styleAxesDaily(ax,xMin,xMax)
% Apply consistent style: daily major grid, LaTeX, Montserrat.

    ax.XLim  = [xMin xMax];
    ax.XTick = xMin:1:xMax;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.GridLineStyle = '-';
    ax.LineWidth = 2.5;
    ax.TickDir = 'out';
    box(ax,'on');
end

function KGE = calcKGE(obs, sim)
% Compute Kling-Gupta Efficiency between obs and sim.
% obs, sim: column vectors (can contain NaNs)

    mask = ~isnan(obs) & ~isnan(sim);
    obs  = obs(mask);
    sim  = sim(mask);

    if numel(obs) < 2
        KGE = NaN;
        return;
    end

    % Components
    r     = corr(obs, sim);                 % linear correlation
    alpha = std(sim) / std(obs);            % variability ratio
    beta  = mean(sim) / mean(obs);          % bias ratio

    KGE = 1 - sqrt( (r - 1).^2 + (alpha - 1).^2 + (beta - 1).^2 );
end
