classdef Initial < handle
  properties
    lb
    ub
    xfunc
    % these are APIs for the Guess() function to pull data out
    r0_guess
    v0_guess
    lat
    incT
  end

  methods
    function x = LatLng_xfunc(obj, r0, v0)
      x = [ r0 v0 ];
    end
  end

  methods(Static)
    function initial = LatLng(launcher, lat, lng);
      body = launcher.body;
      scales = launcher.scales;
      auxdata = launcher.auxdata;
      initial = Initial();
      initial.lat = lat;
      initial.incT = lat;
      initial.r0_guess = auxdata.rbody * [ cos(lat)*cos(lng); cos(lat)*sin(lng); sin(lat) ];
      initial.v0_guess = [ auxdata.omegaMatrix * initial.r0_guess ];
      initial.lb = [ initial.r0_guess' initial.v0_guess' ];
      initial.ub = [ initial.r0_guess' initial.v0_guess' ];
      initial.xfunc = @initial.LatLng_xfunc;
      launcher.initial = initial;
    end
  end
end
