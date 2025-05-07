function [out] = oom(x)
    if x < 1
        o = -ceil(log10(abs(1 ./ x)));
    else
        o = floor(log10(abs(x)));
    end
    out = 10.^o;
end

