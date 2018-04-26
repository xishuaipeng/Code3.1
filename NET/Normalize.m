function data = Normalize(data)
            [row, col] = size(data);
            mean_v = mean(data,1);
            var_v = var(data,[],1);
            mean_v = repmat(mean_v,row,1);
            var_v = repmat(var_v,row,1);
            data = (data - mean_v)./((var_v.^2+ 1e-10).^(0.5));
end
