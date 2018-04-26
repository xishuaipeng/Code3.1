function data = MinMaxData(data)
[row, col] = size(data);
min_v = min(data);
max_v = max(data);
min_v = repmat(min_v,row,1);
max_v = repmat(max_v,row,1);
data = (data - min_v)./(max_v - min_v + 5e-10);
end

