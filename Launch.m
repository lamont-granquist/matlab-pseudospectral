classdef Launch < handle
  properties
    body
    scales
    stages
    auxdata
    terminal
    initial
    result
  end

  methods
    function obj = Launch(body)
      obj.body = body;
    end

    function scales = InitScales(obj, initialMass)
      scales = Scales(obj.body, initialMass);
      obj.scales = scales;
    end

    function auxdata = InitAuxdata(obj, Aref, Cd, Qamax)
      omega = obj.body.omega * obj.scales.time;
      g0 = 9.80665;
      auxdata.omega       = omega;
      auxdata.omegaMatrix = omega * [0 -1 0;1 0 0;0 0 0];
      auxdata.mu          = obj.body.mu/obj.scales.gravparam;
      auxdata.Cd          = Cd;
      auxdata.Aref        = Aref/obj.scales.area;
      auxdata.rho0        = obj.body.rho0 / obj.scales.density;
      auxdata.p0          = obj.body.p0 / obj.scales.pressure;
      auxdata.rbody       = obj.body.radius / obj.scales.length;
      auxdata.H           = obj.body.H / obj.scales.length;
      auxdata.g0          = g0 / obj.scales.accel;
      auxdata.Qamax       = Qamax / obj.scales.pressure; % FIXME: parameter
      auxdata.stages      = obj.stages;
      obj.auxdata = auxdata;
    end

    function dX_dt = VacODE(obj, t, X, thrust, mdot)
      r  = X(1:3);
      v  = X(4:6);
      pv = X(7:9);
      pr = X(10:12);
      m  = X(13);

      u = pv/norm(pv);
      Fm = thrust / m;

      r3 = norm(r)^3;
      r2 = dot(r,r);
      r5 = r2 * r3;

      rdot  = v;
      vdot  = - r / r3 + Fm * u;
      pvdot = - pr;
      prdot = pv / r3 - 3 / r5 * dot(r, pv) * r;

      dX_dt = [ rdot' , vdot', pvdot', prdot', -mdot ]';
    end

    function Guess(obj, problem)
      r0 = obj.initial.r0_guess;
      v0 = obj.initial.v0_guess;
      lat = obj.initial.lat;
      incT = obj.terminal.incT;
      if isstring(incT) && incT == "free_inclination";
        incT = obj.initial.incT;
      end

      % north/east/up
      dr = r0 / norm(r0);
      de = cross([ 0; 0; 1 ], dr);
      de = de / norm(de);
      dn = cross(dr, de);
      dn = dn / norm(dn);

      % launch heading guess from spherical trig
      heading = asin( min(1, max(0, cos(incT)/ cos(lat))) );
      dh = dn * cos(heading) + de * sin(heading);

      % 45 degree pitch up, with a realistic rate vector for a launch
      pitch_guess = 60;
      pv0 = dh * cosd(pitch_guess) + dr * sind(pitch_guess);
      %pr0 = r0 / norm(r0) * 8 / 3;
      pr0 = r0 / norm(r0);

      ti = 0;
      xi = [ r0' v0' pv0' pr0' 0 ];
      for p = 1:length(problem.phases)
        N = problem.phases(p).N;
        nodes = problem.phases(p).nodes;

        if problem.phases(p).auxdata.m0_guess > 0
          xi(13) =  problem.phases(p).auxdata.m0_guess;
        end
        mdot = problem.phases(p).auxdata.mdot;
        thrust = problem.phases(p).auxdata.thrust;
        dt = problem.phases(p).auxdata.burn_guess;
        odefun = @(t,x) obj.VacODE(t, x, thrust, mdot);
        tf = ti + dt;
        range = ( (tf - ti) * nodes + ( tf + ti ) ) / 2;
        [ts, xs] = ode45(odefun, range, xi);
        xi = xs(end,:);

        problem.phases(p).state(1).guess = xs(:,1);
        problem.phases(p).state(2).guess = xs(:,2);
        problem.phases(p).state(3).guess = xs(:,3);
        problem.phases(p).state(4).guess = xs(:,4);
        problem.phases(p).state(5).guess = xs(:,5);
        problem.phases(p).state(6).guess = xs(:,6);
        problem.phases(p).state(7).guess = xs(:,13);
        if thrust > 0
          unorm = vecnorm(xs(:,7:9)')';
          problem.phases(p).control(1).guess = xs(:,7);
          problem.phases(p).control(2).guess = xs(:,8);
          problem.phases(p).control(3).guess = xs(:,9);
          problem.phases(p).control(4).guess = ones(N,1);
        end
        problem.phases(p).guess.time(1) = ti;
        problem.phases(p).guess.time(2) = tf;
        ti = tf;
      end
    end

    function Solve(obj, stages, K, M)
      auxdata = obj.auxdata;
      scales = obj.scales;
      body = obj.body;
      umax = 1 * ones(1,3);
      umin = -1 * ones(1,3);

      p = 1;
      for s = 1:length(stages)
        last_n = stages(s).coasts+1;
        tot = 2 * last_n; % total number of sub-phases in this stage
        if ~stages(s).coast_before
          tot = tot - 1;
        end
        dt = ( stages(s).tf - stages(s).ti ) / tot;
        tstart = stages(s).ti;
        tend = tstart + dt;
        for n = 1:last_n
          if  n > 1 || stages(s).coast_before
            phases(p) = Phase('K', K, 'M', M, 'Nstate', 7, 'Ncontrol', 0, 'Npath', 0);
            phases(p).bounds.initialtime.lower  = -inf;
            phases(p).bounds.initialtime.upper  = inf;
            phases(p).bounds.finaltime.lower    = -inf;
            phases(p).bounds.finaltime.upper    = inf;
            phases(p).bounds.initialstate.lower = [ -inf*ones(1,6) stages(s).mf ];
            phases(p).bounds.initialstate.upper = [ inf*ones(1,6) stages(s).m0 ];
            phases(p).bounds.state.lower        = [ -inf*ones(1,6) stages(s).mf ];
            phases(p).bounds.state.upper        = [ inf*ones(1,6) stages(s).m0 ];
            phases(p).bounds.finalstate.lower   = [ -inf*ones(1,6) stages(s).mf ];
            phases(p).bounds.finalstate.upper   = [ inf*ones(1,6) stages(s).m0 ];
            phases(p).auxdata.mdot = 0;
            phases(p).auxdata.thrust = 0;
            phases(p).auxdata.Ae = 0;
            phases(p).auxdata.jettison = false;
            phases(p).auxdata.stage_num = s;
            phases(p).auxdata.burn_guess = 10 / scales.time;
            phases(p).auxdata.m0_guess = -1;
            p = p + 1;
            tstart = tend;
            tend = tend + dt;
          end

          phases(p) = Phase('K', K, 'M', M, 'Nstate', 7, 'Ncontrol', 4, 'Npath', 1);
          if s == 1 && n == 1
            phases(p).bounds.initialtime.lower  = 0;
            phases(p).bounds.initialtime.upper  = 0;
          else
            phases(p).bounds.initialtime.lower  = -inf;
            phases(p).bounds.initialtime.upper  = inf;
          end
          phases(p).bounds.finaltime.lower    = -inf;
          phases(p).bounds.finaltime.upper    = inf;
          if n == 1
            phases(p).bounds.initialstate.lower = [ -inf*ones(1,6) stages(s).m0 ];
            phases(p).bounds.initialstate.upper = [ inf*ones(1,6) stages(s).m0 ];
          else
            phases(p).bounds.initialstate.lower = [ -inf*ones(1,6) stages(s).mf ];
            phases(p).bounds.initialstate.upper = [ inf*ones(1,6) stages(s).m0 ];
          end
          phases(p).bounds.state.upper        = [ inf*ones(1,6) stages(s).m0 ];
          if n == last_n
            if stages(s).free_final_time
              if stages(s).overburn
                % can burn down the stage to zero
                phases(p).bounds.state.lower        = [ -inf*ones(1,6) 0 ];
                phases(p).bounds.finalstate.lower   = [ -inf*ones(1,6) 0 ];
              else
                % can burn down the stage to mf
                phases(p).bounds.state.lower        = [ -inf*ones(1,6) stages(s).mf ];
                phases(p).bounds.finalstate.lower   = [ -inf*ones(1,6) stages(s).mf ];
              end
              % can't be less than what we start with
              phases(p).bounds.finalstate.upper   = [ inf*ones(1,6) stages(s).m0 ];
            else
              % pinned to exactly end mass
              phases(p).bounds.state.lower        = [ -inf*ones(1,6) stages(s).mf ];
              phases(p).bounds.finalstate.lower   = [ -inf*ones(1,6) stages(s).mf ];
              phases(p).bounds.finalstate.upper   = [ inf*ones(1,6) stages(s).mf ];
            end
          else
            % for interior burns can't overburn the stage, but not pinned to mf
            phases(p).bounds.state.lower        = [ -inf*ones(1,6) stages(s).mf ];
            phases(p).bounds.finalstate.lower   = [ -inf*ones(1,6) stages(s).mf ];
            phases(p).bounds.finalstate.upper   = [ inf*ones(1,6) stages(s).m0 ];
          end
          phases(p).bounds.control.lower      = [ umin 0 ]; % set to 0 to allow throttle down
          phases(p).bounds.control.upper      = [ umax 1 ];
          phases(p).bounds.path.lower         = [ 0 ];
          phases(p).bounds.path.upper         = [ 0 ];
          phases(p).auxdata.mdot = stages(s).mdot;
          phases(p).auxdata.thrust = stages(s).thrust;
          phases(p).auxdata.Ae = stages(s).Ae;
          phases(p).auxdata.stage = stages(s);
          if n == last_n
            phases(p).auxdata.jettison = true;
          else
            phases(p).auxdata.jettison = false;
          end
          phases(p).auxdata.stage_num = s;
          phases(p).auxdata.burn_guess = stages(s).bt / last_n;
          if n == 1
            phases(p).auxdata.m0_guess = stages(s).m0;
          else
            phases(p).auxdata.m0_guess = -1;
          end
          p = p + 1;
          tstart = tend;
          tend = tend + dt;
        end
      end

      e = 1;

      % continuity between phases
      for p = 1:length(phases)-1
        if phases(p).auxdata.jettison
          % splice between a non-coast and the phase after it (jettison)
          eventgroups(e).bounds.lower = [ zeros(1,6) -phases(p).auxdata.stage.mdrop 0 ];
          eventgroups(e).bounds.upper = [ zeros(1,6) -phases(p).auxdata.stage.mdrop 0 ];
        else
          % splice from coast to burn or from a burn without jettison to a coast
          eventgroups(e).bounds.lower = [ zeros(1,6) 0 0 ];
          eventgroups(e).bounds.upper = [ zeros(1,6) 0 0 ];
        end
        e = e + 1;
      end

      % burntime constraints on all phases
      length(phases)
      for p = 1:length(phases)
        if phases(p).auxdata.thrust == 0
          eventgroups(e).bounds.lower = [ 1e-10 ];
          eventgroups(e).bounds.upper = [ inf ];
        else
          if phases(p).auxdata.stage.free_final_time
            eventgroups(e).bounds.lower = [ 1e-10 ];
            if phases(p).auxdata.stage.overburn
              eventgroups(e).bounds.upper = [ phases(p).auxdata.stage.tau * 0.98 ];
            else
              eventgroups(e).bounds.upper = [ phases(p).auxdata.stage.bt ];
            end
          else
            eventgroups(e).bounds.lower = [ phases(p).auxdata.stage.bt ];
            eventgroups(e).bounds.upper = [ phases(p).auxdata.stage.bt ];
          end
        end
        e = e + 1;
      end

      % terminal constraints
      eventgroups(e).bounds.lower = obj.terminal.lb;
      eventgroups(e).bounds.upper = obj.terminal.ub;
      e = e + 1;

      % initial constraint
      eventgroups(e).bounds.lower = obj.initial.lb;
      eventgroups(e).bounds.upper = obj.initial.ub;
      e = e + 1;

      problem.phases      = phases;
      problem.eventgroups = eventgroups;
      problem.auxdata     = auxdata;

      problem.functions.continuous = @obj.launchContinuous;
      problem.functions.endpoint   = @obj.launchEndpoint;
      problem.functions.objective  = @obj.launchObjective;

      problem.auxdata.atm = 1.0;
      problem.auxdata.Qcontrol = 0.0;

      % ipopt options
      opts = struct;
      opts.ipopt.tol = 1e-7; % 1e-7 may cause convergence issues (or may have been poor scaling???)
      opts.ipopt.max_iter = 10000;
      opts.ipopt.print_level = 5;
      opts.ipopt.print_user_options = 'yes';
      opts.ipopt.mumps_permuting_scaling = 7;
      opts.ipopt.mumps_scaling = 8; % this works better than 7 and 77 may not work at all
      opts.ipopt.nlp_scaling_method = 'none';

      solver = Solver();

      solver.precalc(problem);

      obj.Guess(problem);

      solver.gen_guess(problem);

      obj.result = solver.solve(problem, opts);

      opts.lam_x0 = obj.result.lam_x;
      opts.lam_g0 = obj.result.lam_g;

      % warm starting ipopt voodoo: https://list.coin-or.org/pipermail/ipopt/2016-May/004221.html
      %opts.ipopt.warm_start_init_point = 'yes'; % why did I not use this setting from above?
      opts.ipopt.nlp_scaling_method = 'none';
      opts.ipopt.warm_start_bound_frac = 1e-16;
      opts.ipopt.warm_start_bound_push = 1e-16;
      opts.ipopt.warm_start_mult_bound_push = 1e-16;
      opts.ipopt.warm_start_slack_bound_frac = 1e-16;
      opts.ipopt.warm_start_slack_bound_push = 1e-16;

      opts.ipopt.bound_frac = 1e-16; %0.0001;
      opts.ipopt.bound_push = 1e-16; %0.0001;
      opts.ipopt.mu_init = 0.1;
      opts.ipopt.mu_strategy = 'monotone';

      %problem.auxdata.atm = 1.0;

      %solver.copy_guess(obj.result)

      %obj.result = solver.solve(problem, opts);

      %opts.lam_x0 = obj.result.lam_x;
      %opts.lam_g0 = obj.result.lam_g;

      %problem.auxdata.Qcontrol = 1.0;

      %obj.result = solver.solve_guess(obj.result.x, problem, opts);

      obj.plotstuff(obj.result, scales, body);
    end

    function stages = FixStages(obj, stages)
      t = 0;
      for s = 1:length(stages)
        stages(s).ti = t;
        stages(s).tf = t + stages(s).bt;
        t = stages(s).tf;
      end

      for s = 1:length(stages)-1
        stages(s).mdrop = stages(s).mf - stages(s+1).m0;
      end
      stages(end).mdrop = 0;
      obj.stages = stages;
    end

    function plotstuff(obj, result, scales, body)
      result.x(result.problem.phases(end).state(7).idx(end)) * scales.mass
      t0 = result.x(result.problem.phases(1).t0idx) * scales.time;
      tot = result.x(result.problem.phases(end).tfidx) * scales.time - t0
      ct = 0;
      times = [];
      bt = [];
      dv = [];
      for p = 1:length(result.phases)
        dt = (result.x(result.problem.phases(p).tfidx) - result.x(result.problem.phases(p).t0idx)) * scales.time;
        times = [ times dt ];
        if result.phases(p).auxdata.thrust == 0
          ct = ct + dt;
        else
          stage_num = result.phases(p).auxdata.stage_num;
          if stage_num > length(bt)
            bt(stage_num) = 0;
            dv(stage_num) = 0;
          end
          bt(stage_num) = bt(stage_num) + dt;
          tau = result.phases(p).auxdata.stage.tau * scales.time;
          ve = result.phases(p).auxdata.stage.c * scales.speed;
          dv(stage_num) = - ve * log(1 - bt(stage_num) / tau);
        end
      end
      times
      ct
      bt
      sum(bt)
      dv
      sum(dv)

      r = [];
      v = [];
      m = [];
      u = [];
      t = [];
      thr = [];

      for p = 1:length(result.phases)
        r = [ r result.phases(p).state(1:3,:)*scales.length ];
        v = [ v result.phases(p).state(4:6,:)*scales.speed ];
        m = [ m result.phases(p).state(7,:)*scales.mass ];
        if result.phases(p).auxdata.thrust > 0
          u = [ u result.phases(p).control(1:3,:) ];
          thr = [ thr result.phases(p).control(4,:) ];
        else
          N = length(result.phases(p).state(7,:));
          u = [ u zeros(3,N) ];
          thr = [ thr zeros(1,N) ];
        end
        t = [ t result.phases(p).time*scales.time ];
        display_orbit(r(:,end), v(:,end), body);
      end

      m(end) / scales.mass

      O = body.omega * [0 -1 0;1 0 0;0 0 0];
      vsurf = v - O * r;
      h = vecnorm(r) - body.radius;
      rho = body.rho0.*exp(-h/body.H);
      Q = rho.*vecnorm(vsurf).^2/2;

      um = vecnorm(u);
      um(um<=1e-5)=0;
      vsurfm = vecnorm(vsurf);
      cosalpha = dot(u,vsurf)./um./vsurfm;
      cosalpha(isinf(cosalpha)|isnan(cosalpha))=1;
      alpha = acos(cosalpha);
      Qalpha = Q.*alpha;

      figure;
      subplot(3,2,1);
      plot(t,m/1000);
      ylabel('mass, tons', 'fontsize', 14);
      xlabel('Time, sec', 'fontsize', 14);
      hold on;
      grid on;

      subplot(3,2,2);
      plot(t,h/1000);
      ylabel('Height, km', 'fontsize', 14);
      xlabel('Time, sec', 'fontsize', 14);
      hold on;
      grid on;

      subplot(3,2,3);
      plot(t,vecnorm(v)/1000);
      ylabel('velocity, km/s', 'fontsize', 14);
      xlabel('Time, sec', 'fontsize', 14);
      hold on;
      grid on;

      subplot(3,2,4);
      plot(t,thr);
      ylabel('Throttle', 'fontsize', 14);
      xlabel('Time, sec', 'fontsize', 14);
      hold on;
      grid on;

      subplot(3,2,5);
      plot(t,Q);
      ylabel('Dynamic Pressure, Pa', 'fontsize', 14);
      xlabel('Time, sec', 'fontsize', 14);
      hold on;
      grid on;

      subplot(3,2,6);
      plot(t,Qalpha);
      ylabel('Q * alpha, Pa-rads', 'fontsize', 14);
      xlabel('Time, sec', 'fontsize', 14);
      hold on;
      grid on;

    end

    function output = launchEndpoint(obj, input)
      mu = input.auxdata.mu;

      num = 1;

      % continuity constraints between phases
      for i = 1:length(input.phases)-1
        output.eventgroups(num).event = [ input.phases(i+1).initialstate(1:7)-input.phases(i).finalstate(1:7) input.phases(i+1).initialtime-input.phases(i).finaltime ];
        num = num+1;
      end

      % burntime constraints on phases
      for i = 1:length(input.phases)
        output.eventgroups(num).event = [ input.phases(i).finaltime-input.phases(i).initialtime ];
        num = num+1;
      end

      r0 = input.phases(1).initialstate(1:3);
      v0 = input.phases(1).initialstate(4:6);
      rf = input.phases(end).finalstate(1:3);
      vf = input.phases(end).finalstate(4:6);

      % terminal constraint
      output.eventgroups(num).event = obj.terminal.xfunc(rf, vf);
      num = num+1;

      % initial constraints
      output.eventgroups(num).event = obj.initial.xfunc(r0, v0);
      num = num+1;
    end

    function output = launchContinuous(obj, input)
      for p = 1:length(input.phases)
        N     = input.phases(p).N;
        T     = input.phases(p).auxdata.thrust;
        mdot  = input.phases(p).auxdata.mdot;
        Ae    = input.phases(p).auxdata.Ae;
        mu    = input.auxdata.mu;
        O     = input.auxdata.omegaMatrix;
        rbody = input.auxdata.rbody;
        atm   = input.auxdata.atm;
        Qcontrol = input.auxdata.Qcontrol;
        rho0  = input.auxdata.rho0 * atm;
        p0    = input.auxdata.p0 * atm;
        H     = input.auxdata.H;
        Aref  = input.auxdata.Aref;
        Cd    = input.auxdata.Cd;
        Qamax = input.auxdata.Qamax;

        r = input.phases(p).state(:,1:3);
        v = input.phases(p).state(:,4:6);
        m = input.phases(p).state(:,7);
        if T > 0
          u = input.phases(p).control(:,1:3); % pointing direction
          thr = input.phases(p).control(:,4); % throttle [0,1]
        else
          u = zeros(N,3);
          thr = zeros(N,1);
        end


        vsurf = v - r * O';
        vsurfm2 = vsurf(:,1).^2 + vsurf(:,2).^2 + vsurf(:,3).^2 + 1e-8;
        vsurfm = sqrt(vsurfm2);

        rm2 = r(:,1).^2 + r(:,2).^2 + r(:,3).^2 + 1e-8;
        rm = sqrt(rm2);
        rm3 = rm2.^(3/2);

        um2 = u(:,1).^2 + u(:,2).^2 + u(:,3).^2 + 1e-8;
        um = sqrt(um2);

        h = rm - rbody;
        rho = rho0.*exp(-h/H);
        D = - rho*Aref*Cd.*vsurfm.*vsurf/2;
        Tp = p0.*exp(-h/H) * Ae;

%        Q = rho.*vsurfm2/2;

%        udotvs = u(:,1).*vsurf(:,1) + u(:,2).*vsurf(:,2) + u(:,3).*vsurf(:,3);
%
%        % u projected onto vs
%        uvs = udotvs ./ vsurfm2 .* vsurf;
%        uvsm = sqrt( uvs(:,1).^2 + uvs(:,2).^2 + uvs(:,3).^2 );
%        % u perpendicular to vs
%        uvlat = u - uvs;
%        uvlatm = sqrt( uvlat(:,1).^2 + uvlat(:,2).^2 + uvlat(:,3).^2 );
%        alpha = atan2(uvlatm, uvsm);
%
%        cosalpha = udotvs./um./vsurfm;
%        alpha2 = acos(cosalpha);
%        Qalpha = Q.*alpha * Qcontrol / Qamax;

        rdot = v;
        vdot = - mu ./ rm3 .* r + (T*thr-Tp)./m .* u + D./m;
        %mdot = - mdot * ones(size(m)) .* thr;
        mdot = - mdot * ones(size(m));

        output.phases(p).dynamics = [ rdot vdot mdot ];
        if T > 0
          output.phases(p).path = [ um2-1 ];
        end
      end
    end

    function F = launchObjective(obj, input)
      F = 0;
      for p = 1:length(input.phases)
        if input.phases(p).auxdata.jettison % must be last sub-phase of a burn stage
          if input.phases(p).auxdata.stage.free_final_time % must be optimized time
            F = F - input.phases(p).finalstate(7); % optimize for maximum terminal mass
          end
        end
%        r = input.phases(p).state(:,1:3);
%        rm2 = r(:,1).^2 + r(:,2).^2 + r(:,3).^2 + 1e-8;
%        F = F +  sum(penalty_gt(rm2, 1, 1, 1e4));
      end
    end
  end
end
