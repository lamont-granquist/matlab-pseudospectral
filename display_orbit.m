function display_orbit(rf, vf, body)

  [a,eMag,i,O,o,nu,truLon,argLat,lonPer,plr] = rv2orb(rf, vf, body.mu);

  orbit.semi_major_axis = a / 1000;
  orbit.eccentricity = eMag;
  orbit.PeR = (1 - eMag) * a;
  orbit.ApR = (1 + eMag) * a;
  orbit.PeA = orbit.PeR - body.radius;
  orbit.ApA = orbit.ApR - body.radius;
  orbit.inclination = rad2deg(i);
  orbit.LAN = rad2deg(O);
  orbit.argument_of_periapsis = rad2deg(o);
  orbit.true_anomaly = rad2deg(nu);
  orbit.semi_latus_rectum = plr / 1000;

  orbit

end
