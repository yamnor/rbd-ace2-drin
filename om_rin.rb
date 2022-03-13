#!/home/users/sw6/apl/ruby/bin/ruby

#PBS -l select=1:ncpus=1:mpiprocs=1:ompthreads=1:jobtype=core
#PBS -l walltime=2:00:00

pbs_workdir = ENV['PBS_O_WORKDIR']
if pbs_workdir != nil
  Dir.chdir(pbs_workdir)
end

dir = "om_rin"

def cpptraj(cmd)
  File.open("cpp.cmd", "w"){|fw| fw.puts cmd}
  system("module load amber/20/update0; cpptraj < cpp.cmd > cpp.log; rm cpp.cmd cpp.log")
end

def s2n(str)
  return str.split("_").at(1).delete(":").to_i
end

def s2s(str)
  return str.split(":").at(0).downcase.intern
end

def modpdb(pdb_inp, pdb_out)
  pdb = File.open(pdb_inp).readlines
  pdb.map!{|x| x.gsub("CYX", "CYS")}
  File.open(pdb_out, "w"){|fw| fw.puts pdb.join}
end

def exec_ring(dir, nsamples)
  [*1..nsamples].each do |m|
    if File.size("#{dir}/rin_#{m}.pdb") > 100
      modpdb("#{dir}/rin_#{m}.pdb", "#{dir}/rin_#{m}_mod.pdb")
      system("export VICTOR_ROOT=/home/users/sw6/apl/ring/; /home/users/sw6/apl/ring/bin/Ring -i #{dir}/rin_#{m}_mod.pdb --no_energy -E #{dir}/rin_#{m}.edge -N #{dir}/rin_#{m}.node >& #{dir}/rin_#{m}.log")
      system("export VICTOR_ROOT=/home/users/sw6/apl/ring/; /home/users/sw6/apl/ring/bin/Ring -i #{dir}/rin_#{m}_mod.pdb --no_energy --get_iac -t 4.0 -E #{dir}/rin_#{m}.edge_iac -N #{dir}/rin_#{m}.node_iac >& #{dir}/rin_#{m}.log_iac")
    end
  end
end

def read_hbond(filename, nsamples)
  resid, *hbond = File.open(filename).readlines
  resid = resid.split.drop(1).map{|a| x = a.gsub("@", "_").split("_"); *x = [x[1].to_i, x[3].to_i].sort}
  hbond = hbond.map{|a| a.split.drop(1).map(&:to_i)}.transpose
  dat = Array.new(nsamples).map{Array.new}
  nsamples.times do |t|
    resid.zip(hbond).each do |id, hb|
      dat[t] << id if hb[t] > 0
    end
  end
  return dat
end

def read_rin(dir, nsamples)
  vdwrad = {"H": 1.0000, "C": 1.700, "N": 1.625, "O": 1.480, "S": 1.782}
  ritype = {hbond: :hb, vdw: :vdw, ssbond: :ss, ionic: :ion, pipistack: :pp, pication: :pc, iac: :iac}
  dat = Hash.new
  ritype.each_value{|v| dat[v] = Array.new(nsamples).map{Array.new}}
  nsamples.times do |t|
    flg = Array.new(800).map{Array.new(800, true)}
    tmp = File.open("#{dir}/rin_#{t+1}.edge").readlines.drop(1).map{|x| x.split("\t")}.map{|id_1, int, id_2, dist, angle, ene, atom_1, atom_2| [[s2n(id_1), s2n(id_2)].sort, s2s(int), [atom_1[0].intern, atom_2[0].intern], dist.to_f]}
    tmp.each do |id, int, atom, dist|
      if int != :iac
        dat[ritype[int]][t] << id if ritype.has_key?(int)
        flg[id[0]][id[1]] = false
      end
    end
    tmp = File.open("#{dir}/rin_#{t+1}.edge_iac").readlines.drop(1).map{|x| x.split("\t")}.map{|id_1, int, id_2, dist, angle, ene, atom_1, atom_2| [[s2n(id_1), s2n(id_2)].sort, s2s(int), [atom_1[0].intern, atom_2[0].intern], dist.to_f]}
    tmp.each do |id, int, atom, dist|
      if dist - (vdwrad[atom[0]] + vdwrad[atom[1]]) < 0.4
        if flg[id[0]][id[1]]
          dat[ritype[int]][t] << id
        end
      end
    end
  end
  return dat
end

def remove_duplicates(rin, nsamples)
  key = rin.keys
  key.delete(:iac)
  nsamples.times do |t|
    key.each do |k|
      rin[k][t].each do |id|
        rin[:iac][t].delete(id)
      end
    end
  end
  return rin
end

def write_sequence(filename, rin)
  key, val = rin.to_a.transpose
  val = val.transpose
  nsamples = val.size
  File.open(filename, "w") do |fw|
    nsamples.times do |t|
      str = [(t+1).to_s]
      key.each_with_index do |k, i|
        val[t][i].each do |id|
          str << k
          str << id.map(&:to_s).join("\t")
        end
      end
      fw.puts str.join("\t")
    end
  end
end

def calc_fraction(rin)
  key = rin.keys + [:any]
  dim = 800
  num = rin[:hb].size
  tmp = Hash.new
  key.each{|k| tmp[k] = Array.new(num).map{Array.new(dim).map{Array.new(dim, false)}}}
  rin.each do |k, vk|
    vk.each_with_index do |vt, t|
      vt.each do |i, j|
        tmp[k][t][i][j]    = true
        tmp[:any][t][i][j] = true
      end
    end
  end
  dat = Hash.new
  key.each do |k|
    dat[k] = Array.new(dim).map{Array.new(dim, 0.0)}
    for i in 0..(dim-1)
      for j in (i+1)..(dim-1)
        cnt = 0.0
        tmp[k].each do |vt|
          dat[k][i][j] += 1.0 if vt[i][j]
          cnt += 1.0
        end
        dat[k][i][j] /= cnt
      end
    end
  end
  return dat
end

def write_fraction(filename, rin, threshold: 0.05)
  dat = calc_fraction(rin)
  key = dat.keys
  dim = dat[:hb].size
  File.open(filename, "w") do |fw|
    fw.puts ([:i, :j] + key).map(&:to_s).join("\t")
    for i in 0..(dim-1)
      for j in (i+1)..(dim-1)
        if dat[:any][i][j] > threshold
          val = Array.new
          key.each do |k|
            val << dat[k][i][j]
          end
          fw.puts ([i, j] + val.map{|x| "%6.4f"%[x]}).join("\t")
        end
      end
    end
  end
end

cmd = Hash.new

cmd[:hbond] = <<TEXT
parm   om.top
trajin om_1.trj 201 250
trajin om_2.trj 201 250
trajin om_3.trj 201 250
hbond :1-791 angle 117 dist 3.5 avgout #{dir}/hbn.edge series uuseries #{dir}/hbn.series
TEXT

cmd[:pdb] = <<TEXT
parm   om.top
trajin om_1.trj 201 250
trajin om_2.trj 201 250
trajin om_3.trj 201 250
TEXT

nsamples = 150

[*1..nsamples].each do |n|
  cmd[:pdb] += "trajout #{dir}/rin_#{n}.pdb pdb nobox onlyframes #{n}\n"
end

cpptraj(cmd[:hbond])
cpptraj(cmd[:pdb])

exec_ring(dir, nsamples)

rin = read_rin(dir, nsamples)
rin[:hb] = read_hbond("#{dir}/hbn.series", nsamples)

rin = remove_duplicates(rin, nsamples)

write_sequence("#{dir}/rin.sequence", rin)
write_fraction("#{dir}/rin.fraction", rin)
