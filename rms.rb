#!/usr/bin/ruby

cmd = {}

cmd['rmsd'] = <<TEXT
parm [SYS].top
trajin [SYS]_[NUM].trj
rms ToFirst :1-791&!@H= first out [SYS]_[NUM]_rmsd.dat mass
TEXT

cmd['rmsf'] = <<TEXT
parm [SYS].top
trajin [SYS]_[NUM].trj
rms ToFirst :1-791&!@H= first mass
atomicfluct out [SYS]_[NUM]_rmsf.dat byres bfactor
TEXT

cmd['ave'] = <<TEXT
parm [SYS].top
trajin [SYS]_[NUM].trj
rms ToFirst :1-791&!@H= first mass
average [SYS]_[NUM]_ave.pdb pdb
TEXT

['rmsd', 'rmsf', 'ave'].each do |x|
  for sys in ["wt", "om"]
    for num in 1..3
      File.open("cpp.cmd", "w"){|fw| fw.puts cmd[x].gsub("[SYS]", sys).gsub("[NUM]", "%d"%[num])}
      system("/Users/yamnor/.pyenv/shims/cpptraj < cpp.cmd > cpp.log; rm cpp.cmd cpp.log")
    end
  end
end