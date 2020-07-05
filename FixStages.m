function FixStages(stages)
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
end
