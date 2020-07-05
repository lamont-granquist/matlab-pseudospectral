classdef Terminal < handle
  properties
    lb
    ub
    xfunc
    mu
    incT % this is necessary for the initial guess
    LANT % FIXME: implement this for the initial guess for guessing timed launches
  end

  methods
    function x = Kepler5_xfunc(obj, rf, vf)
      hf = cross(rf, vf);
      ef = cross(vf, hf) / obj.mu - rf / norm(rf);
      x = [ hf, ef ];
    end

    function x = FPA3_xfunc(obj, rf, vf)
      sin_gammaf = dot(rf, vf) / norm(rf) / norm(vf);
      x = [ norm(rf) norm(vf) sin_gammaf ];
    end

    function x = FPA4_xfunc(obj, rf, vf)
      sin_gammaf = dot(rf, vf) / norm(rf) / norm(vf);
      hf = cross(rf, vf);
      cosincf = hf(3)/norm(hf);
      x = [ norm(rf) norm(vf) sin_gammaf cosincf];
    end

    function x = FPA5_xfunc(obj, rf, vf)
      sin_gammaf = dot(rf, vf) / norm(rf) / norm(vf);
      hf = cross(rf, vf);
      n = cross([0 0 1], hf);
      cosincf = hf(3)/norm(hf);
      cosLANf = n(1)/norm(n);
      % FIXME: it might be more numerically stable to construct hf here?
      x = [ norm(rf) norm(vf) sin_gammaf cosincf cosLANf ];
    end
  end

  methods(Static)
    function terminal = Kepler5(launcher, smaT, eccT, incT, LANT, ArgPT)
      body = launcher.body;
      scales = launcher.scales;

      terminal = Terminal();
      terminal.mu = body.mu/scales.gravparam;
      [rT, vT] = orb2eci(terminal.mu, [ smaT / scales.length, eccT, incT, ArgPT, LANT, 0 ]);
      hT = cross(rT, vT);
      eT = cross(vT, hT) / terminal.mu - rT / norm(rT);
      terminal.lb = [ hT' eT' ];
      terminal.ub = [ hT' eT' ];
      terminal.incT = incT;
      terminal.LANT = LANT;
      terminal.xfunc = @terminal.Kepler5_xfunc;
      launcher.terminal = terminal;
    end

    function terminal = FPA5(launcher, rT, vT, gammaT, incT, LANT)
      body = launcher.body;
      scales = launcher.scales;

      rT = rT / scales.length;
      vT = vT / scales.speed;

      terminal = Terminal();
      terminal.incT = incT;
      terminal.LANT = LANT;
      terminal.mu = body.mu/scales.gravparam;
      terminal.lb = [ rT vT sin(gammaT) cos(incT) cos(LANT) ];
      terminal.ub = [ rT vT sin(gammaT) cos(incT) cos(LANT) ];
      terminal.xfunc = @terminal.FPA5_xfunc;
      launcher.terminal = terminal;
    end

    function terminal = FPA4(launcher, rT, vT, gammaT, incT)
      body = launcher.body;
      scales = launcher.scales;

      rT = rT / scales.length;
      vT = vT / scales.speed;

      terminal = Terminal();
      terminal.incT = incT;
      terminal.LANT = "free_LAN";
      terminal.mu = body.mu/scales.gravparam;
      terminal.lb = [ rT vT sin(gammaT) cos(incT) ];
      terminal.ub = [ rT vT sin(gammaT) cos(incT) ];
      terminal.xfunc = @terminal.FPA4_xfunc;
      launcher.terminal = terminal;
    end

    function terminal = FPA3(launcher, rT, vT, gammaT)
      body = launcher.body;
      scales = launcher.scales;

      rT = rT / scales.length;
      vT = vT / scales.speed;

      terminal = Terminal();
      terminal.incT = "free_inclination";
      terminal.LANT = "free_LAN";
      terminal.mu = body.mu/scales.gravparam;
      terminal.lb = [ rT vT sin(gammaT) ];
      terminal.ub = [ rT vT sin(gammaT) ];
      terminal.xfunc = @terminal.FPA3_xfunc;
      launcher.terminal = terminal;
    end
  end
end
