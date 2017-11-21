# join-vtk-movie.rb
#
# use as follows:
#
# lfs find . -name '*.*.*.vtk' # check to make sure it's reasonable
#
# mkdir movie
# cp join-vtk-movie.rb join-movie.pbs movie/
# lfs find . -name '*.*.*.vtk' | xargs -I % mv % movie/
#
# cd movie
# sbatch join-movie.pbs
#
system 'mkdir -p merged'

files = Dir.glob('*.P.vtk').map{|f| f.gsub(/.*\.([0-9]{4})\.P.vtk/, '\1')}.uniq

files.each do |num|
  str = Dir.glob("*.#{num}.P.vtk").join(" ")
  cmd = "./join_vtk.x -o merged/P.#{num}.vtk " + str + " >& /dev/null"
  system cmd

  str = Dir.glob("*.#{num}.d.vtk").join(" ")
  cmd = "./join_vtk.x -o merged/d.#{num}.vtk " + str + " >& /dev/null"
  system cmd

  str = Dir.glob("*.#{num}.B1c.vtk").join(" ")
  cmd = "./join_vtk.x -o merged/B1c.#{num}.vtk " + str + " >& /dev/null"
  system cmd

  str = Dir.glob("*.#{num}.B2c.vtk").join(" ")
  cmd = "./join_vtk.x -o merged/B2c.#{num}.vtk " + str + " >& /dev/null"
  system cmd
end
