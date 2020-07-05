classdef Solver < handle
  properties
    guess
    % FIXME: more properties - result, problem, opts, etc.
  end

  methods
    function num = num_constraints(obj, problem);
      num = 0;
      for p = 1:length(problem.phases)
        N = problem.phases(p).N;
        Nstate = problem.phases(p).Nstate;
        Ncontrol = problem.phases(p).Ncontrol;
        Npath = problem.phases(p).Npath;
        num = num + N * ( Nstate + Npath );
      end
      for e = 1:length(problem.eventgroups)
        num = num + length(problem.eventgroups(e).bounds.lower);
      end
    end

    function result = solve(obj, problem, opts);
      lam_x0 = [];
      if isfield(opts, 'lam_x0')
        lam_x0 = opts.lam_x0;
        opts = rmfield(opts, 'lam_x0');
      end

      lam_g0 = [];
      if isfield(opts, 'lam_g0')
        lam_g0 = opts.lam_g0;
        opts = rmfield(opts, 'lam_g0');
      end

      % get the dimensions of the problem
      problem.NConstraints = obj.num_constraints(problem);
      problem.NVariables = length(obj.guess);

      problem.jacobian_row_scales = ones(1,problem.NConstraints);

      % optimizer variables bounds
      [ lbx, ubx ] = obj.gen_bounds(problem);

      % constraint bounds
      [ lbg, ubg ] = obj.gen_constraint_bounds(problem);

      xsym = casadi.SX.sym('x', problem.NVariables);
%      J1 = jacobian(Solver.constraints(xsym, problem), xsym);
%      J = casadi.Function('f', {xsym}, {jacobian(Solver.constraints(xsym, problem), xsym)});
%      J2 = sparse(J(obj.guess));
%      J1(801,1165)
%      obj.guess(725)
%      obj.guess(765)
%      obj.guess(805)
%      obj.guess(1125)
%      obj.guess(965)
%      obj.guess(1085)
%      J2(801,1165)
%      J2(J2 == 0) = 1;
%      max(max(abs(J2)))
%      min(min(abs(J2)))
%      cond(full(J2))


      nlp = struct('x', xsym, 'f', Solver.objective(xsym, problem), 'g', Solver.constraints(xsym, problem));
      solver = casadi.nlpsol('S', 'ipopt', nlp, opts);
      if isempty(lam_x0)
        sol = solver('x0', obj.guess, 'lbx', lbx, 'ubx', ubx, 'lbg', lbg, 'ubg', ubg);
      else
        sol = solver('x0', obj.guess, 'lbx', lbx, 'ubx', ubx, 'lbg', lbg, 'ubg', ubg, 'lam_x0', lam_x0, 'lam_g0', lam_g0);
      end

      result = obj.extract_result(problem, sol);

      %full(Solver.objective(full(sol.x), problem))
      %full(Solver.constraints(full(sol.x), problem))
    end

    function result = extract_result(obj, problem, sol)
      result.problem = problem;
      result.lam_x = sol.lam_x;
      result.lam_g = sol.lam_g;
      x = sol.x;
      result.x = full(x);
      for p = 1:length(problem.phases)
        Nstate = problem.phases(p).Nstate;
        Ncontrol = problem.phases(p).Ncontrol;
        for i = 1:Nstate
          result.phases(p).state(i,:) = full(x(problem.phases(p).state(i).idx));
        end
        for i = 1:Ncontrol
          result.phases(p).control(i,:) = full(x(problem.phases(p).control(i).idx));
        end
        t0 = full(x(problem.phases(p).t0idx));
        tf = full(x(problem.phases(p).tfidx));
        result.phases(p).time = ((tf - t0) * problem.phases(p).nodes + (tf + t0))/2;
        result.phases(p).auxdata = problem.phases(p).auxdata;
      end
    end

    function problem = precalc(obj, problem)
      idx = 1;
      oidx = 1;
      for p = 1:length(problem.phases)
        Nstate = problem.phases(p).Nstate;
        Ncontrol = problem.phases(p).Ncontrol;
        Npath = problem.phases(p).Npath;
        N = problem.phases(p).N;
        for i = 1:Nstate
          last = idx + N-1;
          disp(sprintf("phase %d state %d: %d-%d", p, i, idx, last));
          problem.phases(p).state(i).idx = idx:last;
          idx = last+1;
          olast = oidx + N-1;
          problem.phases(p).dynamics(i).idx = oidx:olast;
          oidx = olast+1;
        end
        for i = 1:Ncontrol
          last = idx + N-1;
          disp(sprintf("phase %d control %d: %d-%d", p, i, idx, last));
          problem.phases(p).control(i).idx = idx:last;
          idx = last+1;
        end
        for i = 1:Npath
          olast = oidx + N-1;
          problem.phases(p).path(i).idx = oidx:olast;
          oidx = olast+1;
        end
        problem.phases(p).t0idx = idx;
        problem.phases(p).tfidx = idx+1;
        idx = idx+2;
      end
    end

    function copy_guess(obj, result)
      guess = result.x;
    end

    function gen_guess(obj, problem)
      for p = 1:length(problem.phases)
        N = problem.phases(p).N;
        Nstate = problem.phases(p).Nstate;
        Ncontrol = problem.phases(p).Ncontrol;

        for i = 1:Nstate
          obj.guess(problem.phases(p).state(i).idx,1) = problem.phases(p).state(i).guess;
        end
        for i = 1:Ncontrol
          obj.guess(problem.phases(p).control(i).idx,1) = problem.phases(p).control(i).guess;
        end
        obj.guess(problem.phases(p).t0idx,1) = problem.phases(p).guess.time(1);
        obj.guess(problem.phases(p).tfidx,1) = problem.phases(p).guess.time(2);
      end
    end

    function [ lb, ub ] = gen_bounds(obj, problem)
      for p = 1:length(problem.phases)
        N = problem.phases(p).N;
        Nstate = problem.phases(p).Nstate;
        Ncontrol = problem.phases(p).Ncontrol;
        for i = 1:Nstate
          lb(problem.phases(p).state(i).idx) = problem.phases(p).bounds.state.lower(i) * ones(N,1);
          ub(problem.phases(p).state(i).idx) = problem.phases(p).bounds.state.upper(i) * ones(N,1);
          if isfield(problem.phases(p).bounds, 'initialstate')
            lb(problem.phases(p).state(i).idx(1)) = problem.phases(p).bounds.initialstate.lower(i);
            ub(problem.phases(p).state(i).idx(1)) = problem.phases(p).bounds.initialstate.upper(i);
          end
          if isfield(problem.phases(p).bounds, 'finalstate')
            lb(problem.phases(p).state(i).idx(end)) = problem.phases(p).bounds.finalstate.lower(i);
            ub(problem.phases(p).state(i).idx(end)) = problem.phases(p).bounds.finalstate.upper(i);
          end
        end
        for i = 1:Ncontrol
          lb(problem.phases(p).control(i).idx) = problem.phases(p).bounds.control.lower(i) * ones(N,1);
          ub(problem.phases(p).control(i).idx) = problem.phases(p).bounds.control.upper(i) * ones(N,1);
          if isfield(problem.phases(p).bounds, 'initialcontrol')
            lb(problem.phases(p).control(i).idx(1)) = problem.phases(p).bounds.initialcontrol.lower(i);
            ub(problem.phases(p).control(i).idx(1)) = problem.phases(p).bounds.initialcontrol.upper(i);
          end
          if isfield(problem.phases(p).bounds, 'finalcontrol')
            lb(problem.phases(p).control(i).idx(end)) = problem.phases(p).bounds.finalcontrol.lower(i);
            ub(problem.phases(p).control(i).idx(end)) = problem.phases(p).bounds.finalcontrol.upper(i);
          end
        end
        lb(problem.phases(p).t0idx) = problem.phases(p).bounds.initialtime.lower;
        lb(problem.phases(p).tfidx) = problem.phases(p).bounds.finaltime.lower;
        ub(problem.phases(p).t0idx) = problem.phases(p).bounds.initialtime.upper;
        ub(problem.phases(p).tfidx) = problem.phases(p).bounds.finaltime.upper;
      end
    end

    function [ cl, cu ] = gen_constraint_bounds(obj, problem)
      cl = [];
      cu = [];
      for e = 1:length(problem.eventgroups)
        cl = [ cl problem.eventgroups(e).bounds.lower ];
        cu = [ cu problem.eventgroups(e).bounds.upper ];
      end
      cl = [ zeros(1,problem.NConstraints-length(cl)) cl ];
      cu = [ zeros(1,problem.NConstraints-length(cu)) cu ];
      for p = 1:length(problem.phases)
        for i = 1:length(problem.phases(p).path)
          cl(problem.phases(p).path(i).idx) = problem.phases(p).bounds.path.lower(i);
          cu(problem.phases(p).path(i).idx) = problem.phases(p).bounds.path.upper(i);
        end
      end
      cl = cl .* problem.jacobian_row_scales;
      cu = cu .* problem.jacobian_row_scales;
    end
  end

  methods(Static)
    function input = marshall_input(x, problem)
      for p = 1:length(problem.phases)
        Nstate = problem.phases(p).Nstate;
        Ncontrol = problem.phases(p).Ncontrol;
        N = problem.phases(p).N;
        if isa(x(1),'casadi.SX')
          input.phases(p).initialstate = casadi.SX.sym('initialstate', 1, Nstate);
          input.phases(p).finalstate = casadi.SX.sym('finalstate', 1, Nstate);
          input.phases(p).state = casadi.SX.sym('state', N, Nstate);
          input.phases(p).control = casadi.SX.sym('control', N, Ncontrol);
        end
        for i = 1:Nstate
          input.phases(p).initialstate(i) = x(problem.phases(p).state(i).idx(1));
          input.phases(p).state(:,i) = x(problem.phases(p).state(i).idx);
          input.phases(p).finalstate(i) = x(problem.phases(p).state(i).idx(end));
        end
        for i = 1:Ncontrol
          input.phases(p).control(:,i) = x(problem.phases(p).control(i).idx);
        end
        input.phases(p).initialtime = x(problem.phases(p).t0idx);
        input.phases(p).finaltime = x(problem.phases(p).tfidx);
        input.phases(p).auxdata = problem.phases(p).auxdata;
        input.phases(p).N = N;
      end

      input.auxdata = problem.auxdata;
    end

    function F = objective(x, problem)
      input = Solver.marshall_input(x, problem);
      F = problem.functions.objective(input);
    end

    function C = constraints(x, problem)
      input = Solver.marshall_input(x, problem);

      C = [];

      output = problem.functions.continuous(input);
      for p = 1:length(output.phases)
        ti = input.phases(p).initialtime;
        tf = input.phases(p).finaltime;
        %disp 'D*state dynamics:'
        %output.phases(p).dynamics(:,7)
        %full(problem.phases(p).D * input.phases(p).state(:,7))
        C = [ C (tf-ti)/2 * output.phases(p).dynamics - problem.phases(p).D * input.phases(p).state output.phases(p).path ];
      end

      C = reshape(C, 1, numel(C));

      output = problem.functions.endpoint(input);
      for e = 1:length(output.eventgroups)
        C = [ C output.eventgroups(e).event ];
      end

      C = C .* problem.jacobian_row_scales;
    end
  end
end
