clear all;
close all;

load('/home/zz111/Spectral_unmixing/cc_diffusion_theory/My_ColorBase_large_shrink_spe_200_200_dt.mat');
data = My_ColorBase_s;
data = full(data);
% %% Step 1: 使用 PCA 将 21 维数据降到 3D
% [coeff, score, ~] = pca(data); % 计算主成分
% reducedData = score(:, 1:3); % 选择前三个主成分
% 
% % 可视化 PCA 结果
% figure;
% scatter3(reducedData(:,1), reducedData(:,2), reducedData(:,3), 5, 'b');
% title('PCA');
% xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
% grid on;
reducedData = data;
%% Method 1
numSamples = 40000;        % size of the 
x = 0.005;                 % threshold of the angle
selectedIdx = [];         % store the idx of the selected points
remainingIdx = 1:size(reducedData, 1);  % The index of remaining points

% ==== 新增初始化部分 ====
minAngles = inf(size(reducedData, 1), 1);  % 初始化最小夹角为无穷大

% Randomly select the initial point
initialIdx = randi(length(remainingIdx), 1);
selectedIdx = [selectedIdx; initialIdx];
remainingIdx(remainingIdx == initialIdx) = [];

for i = 2:numSamples
    i
    if isempty(remainingIdx)
        warning('No more points! Already selected %d points', i-1);
        break;
    end

    % 计算剩余点到最新选中点的夹角
    currentPoint = reducedData(selectedIdx(end), :);
    remainingData = reducedData(remainingIdx, :);
    
    % 计算点积和模长
    dotProduct = remainingData * currentPoint';
    normCurrent = norm(currentPoint);
    normRemaining = vecnorm(remainingData, 2, 2);
    
    % 数值稳定性处理
    epsilon = 1e-10;
    cosineSim = dotProduct ./ (normRemaining * normCurrent + epsilon);
    cosineSim = min(max(cosineSim, -1), 1);  % 限制在 [-1, 1] 范围内
    newAngles = acos(cosineSim);             % 计算夹角（弧度）

    minAngles(remainingIdx) = min(minAngles(remainingIdx), newAngles);
    
    % 筛选最小夹角 > x 的候选点
    validCandidates = remainingIdx(minAngles(remainingIdx) > x);
    
    if isempty(validCandidates)
        warning('无更多满足夹角阈值 x=%.2f 的点！已选中 %d 个点', x, i-1);
        break;
    end
    
    % 选择夹角最大的候选点
    [~, farthestIdx] = max(minAngles(validCandidates));
    nextIdx = validCandidates(farthestIdx);
    
    % Update the selected index and remaing index
    selectedIdx = [selectedIdx; nextIdx];
    remainingIdx(remainingIdx == nextIdx) = [];
end

% Extract the subset
subsetData = data(selectedIdx, :);

% Visulization
figure;
scatter3(reducedData(:,1), reducedData(:,2), reducedData(:,3), 5, 'b'); hold on;
scatter3(reducedData(selectedIdx,1), reducedData(selectedIdx,2), reducedData(selectedIdx,3), 20, 'r', 'filled');
title(sprintf('PCA + 距离阈值采样 (x=%.2f)', x));
legend('所有数据', '选中的子集');
grid on;


%% 选择使用哪种方法的结果
save('subset_cita_distance0.005.mat', 'subsetData');
