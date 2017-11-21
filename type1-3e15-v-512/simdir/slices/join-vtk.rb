# join-vtk.rb: calls join_vtk.x where needed
#
# This script joins vtk files as a simulation runs.  Assumes the
# raw vtk files live in directories id*/*.vtk and puts the merged
# files in merged/*.vtk.
#
# Only merges files which don't exist or are out of date, so can be
# run progressively as a simulation runs.
#
# Note that this does *not* work with SMR.
#
require 'fileutils'

def issue_cmd(cmd)
  system cmd
end

def strip_digits(str)
  str.gsub(/.*\.([0-9]{4})\..*/, '\1')
end

def get_base(str)
  front = str.gsub(/(.*)\.([0-9]{4})\..*/, '\1') # eg, .dddd.vtk
  front.sub(/.*\//, '')                          # eg, id15/lev1/
end


if not File.executable?('./join_vtk.x')
  puts "###error: join_vtk.x not found"
  exit 1
end


issue_cmd 'mkdir -p merged'

dirs  = Dir.glob('id[^0]*').map{|f| f.sub(/^id/, '')}
files = Dir.glob('id0/*.vtk').map{|f| strip_digits(f)}
base  = get_base(Dir.glob('id0/*.vtk').first)

files.each do |num|
  outfile = "merged/#{base}.#{num}.vtk"

  infiles = ["id0/#{base}.#{num}.vtk"]
  infiles = infiles + dirs.map{|d| "id#{d}/#{base}-id#{d}.#{num}.vtk"}

  unless FileUtils.uptodate?(outfile, infiles)
    puts "writing #{outfile}..."
    str = infiles.join(" ")
    cmd = "./join_vtk.x -o #{outfile} " + str + '>& /dev/null'
    issue_cmd cmd
  end
end


file1 = "id0/#{base}.hst"
file2 = "merged/#{base}.hst"

unless FileUtils.uptodate?(file1, file2)
  puts "copying hst file..."
  cmd = "cp #{file1} #{file2}"
  issue_cmd cmd
end

exit 0