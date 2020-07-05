classdef Scales < handle
  properties(SetAccess=immutable)
    length
    speed
    time
    accel
    mass
    force
    mdot
    area
    volume
    density
    pressure
    gravparam
  end

  methods
    function obj = Scales(body, m0)
      obj.length    = body.radius;
      obj.speed     = sqrt( body.mu / obj.length );
      obj.time      = obj.length / obj.speed;
      obj.accel     = obj.speed / obj.time;
      obj.mass      = m0;
      obj.force     = obj.mass * obj.accel;
      obj.mdot      = obj.mass / obj.time;
      obj.area      = obj.length^2;
      obj.volume    = obj.area * obj.length;
      obj.density   = obj.mass / obj.volume;
      obj.pressure  = obj.force / obj.area;
      obj.gravparam = obj.accel * obj.length^2;
    end
  end
end
