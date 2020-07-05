classdef Stage < handle
  properties
    m0
    mf
    bt
    thrust
    mdot
    c
    a0
    tau
    Ae
    isp_vac
    isp_atm
    coasts
    coast_before
    ti
    tf
    mdrop
    free_final_time
    overburn
    scales
  end

  methods
    function obj = Stage(scales, thrust, m0, mf, bt)
      obj.m0 = m0 / scales.mass;
      obj.mf = mf / scales.mass;
      obj.bt = bt / scales.time;
      obj.thrust = thrust / scales.force;
      obj.mdot = ( ( m0 - mf ) / bt ) / scales.mdot;
      obj.c = obj.thrust / obj.mdot; % normalized ve
      obj.a0 = obj.thrust / obj.m0; % starting accel
      obj.tau = obj.m0 / obj.mdot; % burn to zero time
      obj.coast_before = false;
      obj.coasts = 0;
      obj.Ae = 0;
      obj.free_final_time = false;
      obj.scales = scales;
      obj.overburn = false;
    end

    function SetISP(obj, isp_vac, isp_atm)
      obj.isp_vac = isp_vac;
      obj.isp_atm = isp_atm;
      obj.Ae = exit_area(obj.mdot * obj.scales.mdot, isp_vac, isp_atm) / obj.scales.area;
    end

    function SetAe(obj, Ae)
      obj.Ae = Ae / obj.scales.area;
    end
  end

  methods(Static)
    function stage = StageFromMdot(scales, thrust, m0, mdot, bt)
      mf = m0 - mdot * bt;
      stage = Stage(scales, thrust, m0, mf, bt);
    end

    function stage = StageFromKSP2(scales, twr, m0, mf, isp_vac, isp_atm)
      g0 = 9.80665;
      thrust = twr * m0 * g0
      mdot = thrust / (isp_vac * g0 );
      bt = (m0 - mf) / mdot;
      stage = Stage(scales, thrust, m0, mf, bt);
      stage.SetISP(isp_vac, isp_atm);
    end

    function stage = StageFromKSP(scales, twr, slt, m0, mf, bt)
      g0 = 9.80665;
      mdot = (m0 - mf) / bt;
      thrust = twr * m0 * g0
      isp_vac = thrust / mdot / g0;
      isp_atm = isp_vac * slt/twr;
      stage = Stage(scales, thrust, m0, mf, bt);
      stage.SetISP(isp_vac, isp_atm);
    end
  end
end


