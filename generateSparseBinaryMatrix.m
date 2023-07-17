function binaryMatrix = generateSparseBinaryMatrix(rows, cols, sparsity)
    % 计算矩阵中应该有多少非零元素
    numNonZeros = round(rows * cols * sparsity);

    % 生成随机的非零元素的索引
    indices = randperm(rows * cols, numNonZeros);

    % 创建一个全零矩阵
    binaryMatrix = zeros(rows, cols);

    % 将随机选择的索引位置设置为 1
    binaryMatrix(indices) = 1;
end
