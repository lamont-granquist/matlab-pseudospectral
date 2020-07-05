classdef Body < handle
  properties
    radius
    period
    mu
    omega
    rho0
    p0
    H
  end

  methods(Static)
    function body = Earth()
      body = Body();
      body.radius = 6378145;
      body.period = 86164.091;
      body.mu     = 3.986012e14;
      body.omega  = 7.29211585e-5;
      body.rho0   = 1.225;
      body.p0     = 101325;
      body.H      = 7200;
    end

    function body = Kerbin()
      body = Body();
      body.radius = 600000;
      body.period = 21549.425;
      body.mu     = 3.5316000e12;
      body.omega  = 2.91570904e-4;
      body.rho0   = 1.225;
      body.p0     = 101325;
      body.H      = 5760;
    end
  end
end


