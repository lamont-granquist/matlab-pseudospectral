classdef Phase < handle
  properties(SetAccess=immutable)
    K
    M
    D
    N
    weights
    nodes
    Nstate
    Ncontrol
    Npath
  end

  properties
    state
    control
    bounds
    dynamics
    path
    guess
    tfidx
    t0idx
    auxdata
  end

  methods
    function obj = Phase(varargin)
      p = inputParser;
      p.addParameter('K', 0, @(x) (x == floor(x) && x > 0));
      p.addParameter('M', 0, @(x) (x == floor(x) && x > 0));
      p.addParameter('Nstate', 0, @(x) (x == floor(x) && x >= 0));
      p.addParameter('Ncontrol', 0, @(x) (x == floor(x) && x >= 0));
      p.addParameter('Npath', 0, @(x) (x == floor(x) && x >= 0));
      p.parse(varargin{:});
      fn = fieldnames(p.Results);
      for k = 1:numel(fn)
        eval(sprintf('obj.%s = %d;', fn{k}, p.Results.(fn{k})));
      end

      %[ obj.D, obj.nodes, obj.weights, obj.N ] = d_matrix2(obj.M, obj.K, 2);
      [ obj.D, obj.nodes, obj.weights, obj.N ] = d_matrix(obj.M, obj.K);

      state = Variable.empty
      for i = 1:obj.Nstate
        state(i) = Variable(obj.N);
      end
      control = Variable.empty
      for i = 1:obj.Ncontrol
        control(i) = Variable(obj.N);
      end
      obj.state = state;
      obj.control = control;
    end

    function guess_api(obj, what, val)
      switch what
        case 'constant_control'
          if length(val) ~= obj.Ncontrol
            error( "bad length" );
          end
          for i = 1:obj.Ncontrol
            obj.control(i).constant_guess(val(i))
          end
        case 'constant_state'
          if length(val) ~= obj.Nstate
            error( "bad length" );
          end
          for i = 1:obj.Nstate
            obj.state(i).constant_guess(val(i))
          end
        case 'time'
          if length(val) ~= 2
            error( "bad length" );
          end
          obj.guess.time = val';
        otherwise
          error( sprintf("unknown option %s", what) );
      end
    end
  end

  methods(Static)
    function phases = PhaseArray(n, varargin)
      for i = 1:n
        phases(i) = Phase(varargin{:});
      end
    end

    function phases = PhasesFromStages(stages)
    end
  end
end
