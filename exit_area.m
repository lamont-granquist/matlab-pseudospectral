function Ae = exit_area(mdot, ispvac, ispatm)
  g0 = 9.80665;
  p0 = 101325;
  Ae = g0 / p0 * mdot * ( ispvac - ispatm);
end
