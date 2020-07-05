function F = penalty_gt(x, a, penalty, sharpness)
   F = sigmoid((x - a) .* sharpness) .* ((x - a) + penalty);
end
