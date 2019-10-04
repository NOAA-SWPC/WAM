#!/usr/bin/ruby

RSYNC="/usr/bin/rsync -vrLpgoDu -e ssh"
#BASEDIR="/lfs0/projects/rtfim/FIM/FIMrun"
BASEDIR="/lfs0/projects/rtfim/FIMX/FIMrun"
DEBUG=true
KEEP=2

# Get most recent run
now=Time.now.to_i
last_run=Time.at(now - now % (3600*12))

0.upto(KEEP-1) { |i|

  yyyymmddhhmm=(last_run - (i*3600*12)).strftime("%Y%m%d%H%M")
  hh=(last_run - (i*3600*12)).strftime("%H").to_i

  Dir["#{BASEDIR}/fim_[0-9]*_[0-9]*_[0-9]*_#{yyyymmddhhmm}/post_C/fim/NAT/grib1"].each { |dir|
    cmd="#{RSYNC} #{dir}/[0-9]* /rt/fimx/nat/grib1"
    output=`#{cmd} 2>&1`
    error=$?.exitstatus
    puts output if DEBUG
    if error != 0
      puts "ERROR: #{cmd} failed!  Exit status=#{error}"
      puts output
    end

  }

  # Only rsync grib2 files that are valid at 00Z and 12Z
  Dir["#{BASEDIR}/fim_[0-9]*_[0-9]*_[0-9]*_#{yyyymmddhhmm}/post_C/fim/NAT/grib2"].each { |dir|
    files12Z=Dir["#{dir}/[0-9]*[0-9]"].delete_if {|file| (file.slice(-4,4).to_i % 12) != 0 }.join(" ")
    cmd="#{RSYNC} #{files12Z} /rt/fimx/nat/grib2"
    output=`#{cmd} 2>&1`
    error=$?.exitstatus
    puts output if DEBUG
    if error != 0
      puts "ERROR: #{cmd} failed!  Exit status=#{error}"
      puts output
      next
    end

  }

  if File.exists?("#{BASEDIR}/fim_7_64_240_#{yyyymmddhhmm}/tracker_C/168")
    cmd="#{RSYNC} --include '*track*' --include '*.grib.F7C' --include '*/' --exclude '*' #{BASEDIR}/fim_7_64_*_#{yyyymmddhhmm}/tracker_C/168/* /rt/fimx/tracker"
    output=`#{cmd} 2>&1`
    error=$?.exitstatus
    puts output if DEBUG
    if error != 0
      puts "ERROR: #{cmd} failed!  Exit status=#{error}"
      puts output
    end
  end

}
