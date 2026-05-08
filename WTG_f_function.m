function [f] = WTG_f_function(x)

f = x.* asinh(x) - sqrt(1 + x.^2);

end

