clear;
clc;
close all;

%% Model and numerical parameters
nRows = 5;
nCols = 5;

epsilon = 0.02;
f       = 1.4;
q       = 0.002;
Du      = 0.20;

dx      = 1.0;
dt      = 1.0e-3;
tFinal  = 10.0;        % <-- CHANGE THIS to see more/less time
nSteps  = floor(tFinal / dt) + 1;
time    = (0:nSteps-1) * dt;

plotEvery = 200;

matFileName = 'oregonator_5x5_results.mat';
csvFileName  = 'oregonator_5x5_results.csv';

%% Initial conditions
rng(1);
u = 0.20 + 0.02 * randn(nRows, nCols);
v = 0.10 + 0.02 * randn(nRows, nCols);

u(3,3) = 0.80;   % stronger initial kick so red appears sooner
u = max(u, 1e-6);
v = max(v, 1e-6);

%% Preallocate storage
U_history = zeros(nRows, nCols, nSteps);
V_history = zeros(nRows, nCols, nSteps);
U_history(:,:,1) = u;
V_history(:,:,1) = v;

%% Live heatmap figure
figure('Name','Oregonator 5x5 Heatmap','Color','w');
hImg = imagesc(u);
axis equal tight;
colormap(turbo);
colorbar;
caxis([0 1]);
title(sprintf('Activator Concentration u at t = %.3f', 0));
xlabel('Column Index');
ylabel('Row Index');
set(gca,'XTick',1:nCols,'YTick',1:nRows,'LineWidth',1.0,...
    'XGrid','on','YGrid','on','GridColor','k','GridAlpha',0.9);
hold on;
for xLine = 0.5:1:(nCols+0.5)
    plot([xLine xLine],[0.5 nRows+0.5],'k-','LineWidth',1.1);
end
for yLine = 0.5:1:(nRows+0.5)
    plot([0.5 nCols+0.5],[yLine yLine],'k-','LineWidth',1.1);
end
hold off;
drawnow;

%% Time integration (Forward Euler)
for step = 2:nSteps
    lapU = neumannLaplacian(u, dx);   % calls the separate file

    reactionU = (1/epsilon) .* (u - u.^2 - f.*v.*((u - q)./(u + q)));
    reactionV = u - v;

    u = max(u + dt*(reactionU + Du*lapU), 0);
    v = max(v + dt*reactionV,             0);

    U_history(:,:,step) = u;
    V_history(:,:,step) = v;

    if mod(step, plotEvery)==0 || step==nSteps
        set(hImg,'CData',u);
        title(sprintf('Activator Concentration u at t = %.3f', time(step)));
        drawnow;
    end
end

%% Time-series analysis for center cell (3,3)
centerRow = 3; centerCol = 3;
uCenter = squeeze(U_history(centerRow,centerCol,:));
vCenter = squeeze(V_history(centerRow,centerCol,:));

transientCut   = round(0.20*nSteps);
analysisTime   = time(transientCut:end);
analysisSignal = uCenter(transientCut:end);

peakMask = false(size(analysisSignal));
for k = 2:numel(analysisSignal)-1
    if analysisSignal(k) > analysisSignal(k-1) && ...
       analysisSignal(k) >= analysisSignal(k+1) && ...
       analysisSignal(k) > mean(analysisSignal)
        peakMask(k) = true;
    end
end

peakIndicesLocal = find(peakMask);
peakTimes = analysisTime(peakIndicesLocal);

if numel(peakTimes) >= 2
    meanPeriod           = mean(diff(peakTimes));
    oscillationFrequency = 1 / meanPeriod;
else
    meanPeriod = NaN; oscillationFrequency = NaN;
end

fprintf('Center cell (%d,%d) frequency = %.6f\n', centerRow, centerCol, oscillationFrequency);
fprintf('Estimated period              = %.6f\n', meanPeriod);
fprintf('Detected peaks after transient= %d\n',  numel(peakTimes));

%% Time-series plot
figure('Name','Time Series for Cell (3,3)','Color','w');
plot(time, uCenter,'b-','LineWidth',1.5); hold on;
plot(time, vCenter,'r--','LineWidth',1.2);
if ~isempty(peakTimes)
    plot(peakTimes, analysisSignal(peakIndicesLocal),...
         'ko','MarkerFaceColor','y','MarkerSize',5);
end
hold off; grid on;
xlabel('Time'); ylabel('Concentration');
title(sprintf('Cell (3,3): Freq = %.6f', oscillationFrequency));
legend('u(3,3)','v(3,3)','Peaks','Location','best');

%% Save results
results.time                 = time;
results.U_history            = U_history;
results.V_history            = V_history;
results.u_center             = uCenter;
results.v_center             = vCenter;
results.peak_times           = peakTimes;
results.oscillation_frequency= oscillationFrequency;
results.mean_period          = meanPeriod;
results.parameters = struct('epsilon',epsilon,'f',f,'q',q,'Du',Du,...
    'dx',dx,'dt',dt,'tFinal',tFinal,'nRows',nRows,'nCols',nCols);
save(matFileName,'results');
fprintf('Saved: %s\n', matFileName);
