classdef Variable < handle
  properties(SetAccess=immutable)
    N
  end

  properties
    guess
    idx
    bounds
  end

  methods
    function obj = Variable(N)
      obj.N = N;
    end

    function constant_guess(obj, val)
      obj.guess = val * ones(obj.N, 1);
    end
  end
end
