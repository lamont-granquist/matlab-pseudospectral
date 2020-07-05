function F = penalty_gt(x, a, penalty, sharpness)
   F = sigmoid((a - x) .* sharpness) .* ((a - x) + penalty);
end
