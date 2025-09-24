function showProgress(iter, total)
% showProgress 显示命令行进度条（兼容 Live Script / Command Window）
%
%   showProgress(i, N)
%   i     - 当前迭代次数
%   N     - 总迭代次数
%
% 示例:
%   N = 100;
%   for i = 1:N
%       pause(0.05); % 模拟计算
%       showProgress(i, N);
%   end

    % 进度百分比
    percent = iter / total;
    barLen  = 50;                           % 进度条长度
    numBars = round(barLen * percent);
    barStr  = sprintf('[%s%s] %3.0f%% (%d/%d)', ...
                      repmat('=',1,numBars), repmat(' ',1,barLen-numBars), ...
                      percent*100, iter, total);

    % 使用 persistent 保存上一次输出（解决 Live Script 不识别 \r 的问题）
    persistent lastStr;
    if isempty(lastStr)
        lastStr = '';
    end

    % 删除上一行（适配 Live Script）
    fprintf(repmat('\b',1,length(lastStr)));
    fprintf('%s', barStr);

    % 保存当前字符串，供下次覆盖
    lastStr = barStr;

    % 循环结束时换行并清空缓存
    if iter == total
        fprintf('\n');
        lastStr = '';
    end
end