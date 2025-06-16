function [] = plot_initial_ensemble(datasets,xlabels,nbins)

% Example datasets (replace these with your actual data)
% dataset1 = randn(1000,1);  % Normal distribution (Example data)
% dataset2 = randn(1000,1) * 2 + 5;  % Another normal distribution
% dataset3 = exprnd(1, 1000, 1);  % Exponential distribution (Example data)


% Number of datasets
n = size(datasets,2);

% Calculate the number of rows and columns for subplots
if n < 3
    nrows = 1;
    ncols = n;  % 1 row for less than 3 datasets
else
    nrows = ceil(n / 3);  % Multiple rows if n >= 3
    ncols = 3;  % Fix number of columns to 3
end

% Create a figure
figure;

% Loop over datasets to create subplots
for i = 1:n
    % Create a subplot for the i-th dataset
    ax1 = subplot(nrows, ncols, i);
    hold on;
    
% Number of datasets
n = size(datasets,2);

% Calculate the number of rows and columns for subplots
if n < 3
    nrows = 1;
    ncols = n;  % 1 row for less than 3 datasets
else
    nrows = ceil(n / 3);  % Multiple rows if n >= 3
    ncols = 3;  % Fix number of columns to 3
end

% Loop over datasets to create subplots
for i = 1:n
    % Create a subplot for the i-th dataset
    ax = subplot(nrows, ncols, i);
    hold on;
    
    % Get current dataset
    dataset = datasets(:,i);
    
        % Calculate number of bins automatically using the 'auto' method
    numBins = nbins;  % You can adjust this or set it dynamically if needed
    
    % Plot the histogram with proper normalization and appearance
    h = histogram(dataset, numBins, 'FaceColor', [0.2, 0.4, 0.6], 'EdgeColor', 'none');
    
    % Compute the sample mean and standard deviation
    sampleMean = mean(dataset);
    sampleStd = std(dataset);
    
    % Add labels to the plot
    if i == 1
         text(mean(dataset), max(h.Values) * 0.8, ...
         sprintf('Bins: %d\nMean: %.7f\nStd: %.7f', numBins, sampleMean, sampleStd), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
         'FontSize', 12, 'FontName', 'Garamond', 'Interpreter', 'latex'); 
         % xscale log
    else

    
    text(mean(dataset), max(h.Values) * 0.8, ...
         sprintf('Bins: %d\nMean: %.2f\nStd: %.2f', numBins, sampleMean, sampleStd), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
         'FontSize', 12, 'FontName', 'Helvetica', 'Interpreter', 'latex');   
    if i == 6
        % xscale log
    end
    end
    % Set axis labels with LaTeX formatting
    xlabel(ax, xlabels{i}, 'Interpreter', 'latex', 'FontSize', 16);
    ylabel(ax, 'Count', 'Interpreter', 'latex', 'FontSize', 16);
    
    % Set font size for tick labels
    set(ax, 'FontSize', 14, 'FontName', 'Garamond');

    % 
    if i == 1
        % xscale log
    end
    
    % Set tick direction to outside
    set(ax, 'TickDir', 'out');
    
    % Title (can also include the number of bins)
    % title(ax, sprintf('Dataset %d', i), 'Interpreter', 'latex', 'FontSize', 16);
    
    % Set grid
    grid on;
    hold off;
end

% Adjust figure layout to fit within A4 paper size
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [21, 29.7]);  % A4 size paper
set(gcf, 'PaperPosition', [0, 0, 21, 14]);  % Fit to A4 width

end

