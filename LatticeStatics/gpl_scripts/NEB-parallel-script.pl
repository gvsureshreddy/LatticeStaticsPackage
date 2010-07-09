#! /usr/bin/perl
#use Term::ReadKey;
use POSIX ":sys_wait_h";
use Cwd;
use FileHandle;
use File::Basename;
use File::Copy;
use File::Find;

# argument list
#  - wall-clock time (when should I START shutting down)
#  - relative path to executable (../column)
#  - relative path to root directory of bfb-tree

$starttime = time();

if ((scalar @ARGV) == 0)
{
  die "Usage\n" .
        "NEB-auto <wall-clock time in seconds for when to START shutting down> <rel-path-name-of-executable> <rel-path-to-NEBrootdir>\n";
}

$walltimeout = shift;
$ProgExec = getcwd() . "/" . shift;
$RootBFBDir = getcwd() . "/" . shift;

$PROCESSFLAG = "WAITING";
$ABORTEDFLAG = "ABORTED";

STDIN->autoflush(1);
print "\n\nStarting Parallel QCNEB\n\n";

# global variables
$geo = "columnNi.geo";
$pots = "Potentials";
$sleeptime = 15;
#
$endtime = $starttime + $walltimeout;
$decrement = ceiling($walltimeout/100);

############ get processor names #################
$i = 0;
$TotalNumProcessors = 0;
@ProcList = glob "/tmp/*.q/machines";
if (scalar @ProcList != 1)
{
  die "Could not find unique 'machines' file.\n;";
}
else
{
  print "Getting processor list from file: $ProcList[0]\n";
}

open(LST,"<$ProcList[0]");
while(<LST>)
{
  if (!/^$/)
  {
    chomp $_;
    if ($i > 0)
    {
      push @cpulist, $_;
      print "Processor $i has name: ",$cpulist[-1],"\n";
      $TotalNumProcessors++;
    }
    else
    {
      print "Processor $i has name: ",$_, " and is being reserved for this script's use.\n",
            "\t\t(it is assumed that I am already running on this processor....)\n";
    }
    $i++;
  }
}
close(LST);
print "Total number of processors available for use is $TotalNumProcessors.\n\n";
############ create runtime sentinel and timer ######################

$maintimerfile = "$RootBFBDir/#MASTER-TIMER#";
print "Setting wall-clock timeout to ", scalar localtime($endtime), ".... ";
defined($TimerPID = fork()) or die "can't fork timer.\n";
if ($TimerPID == 0)
{
  # Timer's job
  $remainingtime = $endtime - time();
  if ($decrement > $remainingtime)
  {
    $decrement = $remainingtime;
  }
  open(MAINTIMEFILE,">$maintimerfile");
  MAINTIMEFILE->autoflush(1);

  print MAINTIMEFILE "$remainingtime seconds (",$remainingtime/60," minutes) remaining at ", scalar localtime(time()), ".\n";
  while(time() < $endtime)
  {
    sleep($decrement);
    if (! -e $maintimerfile) # check if someone removed the timer file
    {
      exit;
    }
    $remainingtime = $endtime - time();
    print MAINTIMEFILE "$remainingtime seconds (",$remainingtime/60," minutes) remaining at ", scalar localtime(time()), ".\n";
  }

  close(MAINTIMEFILE);

  unlink("$maintimerfile");

  exit;
}
print "Done. Timer pid = $TimerPID.\n\n";


########### clean up aborted directories and files ###############
print "Processing $ABORTEDFLAG directories.\n";
$i=0;


@abtd = ();
find_names($RootBFBDir,"#$ABORTEDFLAG#",\@abtd);
foreach (@abtd)
{
  $i++;
  $workingdir = $_;
  $workingdir =~ s/\/#$ABORTEDFLAG#//;

  @dummy = split('/',$workingdir); 
  $shortdir = pop @dummy;

  print "     Working on: ", $shortdir,"\n";
  unlink glob("$workingdir/*.plt.gz");
  unlink "$workingdir/abort.dat";
  unlink glob("$workingdir/*.out.gz");
  unlink "$workingdir/qc.{log,cmd}.gz";

  move("$workingdir/#$ABORTEDFLAG#","$workingdir/#$PROCESSFLAG#");
}
if ($i == 0)
{
  print "     No $ABORTEDFLAG directories found.\n";
}
print "Done processing $ABORTEDFLAG direcotries.\n\n";


####################### start running #####################
@runningjobs = ();
while((-e $maintimerfile) &&
      ( ((scalar @cpulist) < $TotalNumProcessors) ||
        (defined(find_first($RootBFBDir,"#$PROCESSFLAG#")))
      )
     )
{
  #### find any jobs that have finished from the 
  #### list of running job directories (@pths)
  clean_up_finished_jobs();

  #### Start as many jobs as possible
  start_waiting_jobs();

  #### if no more waiting jobs are available sleep
  if (scalar @cpulist > 0)
  {
    # there are no paths ready to be computed then sleep
    sleep($sleeptime);
  }
}
# finished main loop

# set abort.dat for any running processes
foreach(@runningjobs)
{
  $abortdir = $_;

  $flnm = $abortdir;
  $flnm =~ s/.*\/([^\/]*)$/$1/;
  print "Aborting $flnm.\n";
  open(ABORTFL,">$abortdir/abort.dat");
  print ABORTFL "T\n";
  close(ABORTFL);
}
print "\n";

# check if time remains
if (-e "$maintimerfile")
{
  print "Finished and Timer file still exists.  Removing.\n";
  kill 9, $TimerPID;
  unlink "$maintimerfile";
}

# wait for all child processes to exit
print "Waiting for child processes to end...\n";
#### find any processors that have finished
while ( (scalar @runningjobs) > 0)
{
  clean_up_finished_jobs();

  sleep($sleeptime);
}

while (-1 != ($ans =wait))
{
  print "     Process $ans ended.\n";
}
print "Done.\n\n";

# clean up

print "Ending Parallel QCNEB.\n\n";
exit;

#--------------------------------------------------------------------------
sub clean_up_finished_jobs
{
  # uses globals:
  #   runningjobs
  #   cpulist
  #   RunningProcesses
  #
  my @remainingjobs = ();
  foreach (@runningjobs)
  {
    my $curdir = $_;
    my @tmp = split('/',$curdir);
    my $flnm = pop @tmp;

    if (! -e "$curdir/#PROCESSING#")
    {
      #### job has finished/exited
      # clean up after the run returns.
      
      # use -q to stop warnings
      system("gzip -f $curdir/$flnm.in $curdir/$flnm.out $curdir/*.plt" .
         " $curdir/qc.{log,cmd} >& /dev/null");

      # add cpu for current job to available cpulist
      push @cpulist, $RunningProcesses{$curdir};
      # remove job key from runningprocesses hash
      delete $RunningProcesses{$curdir};
      # print message to screen based on exit status of job
      if (-e "$curdir/#ABORTED#")
      {
        print "Process aborted on $RunningProcesses{$curdir} at ",scalar localtime(time()),
        " (", (scalar @cpulist)," processors available).\n\n";
      }
      elsif (-e "$curdir/#ERROR#")
      {
        print "Process exited with ERROR on $RunningProcesses{$curdir} at ",
        scalar localtime(time()), " (", (scalar @cpulist)," processors available).\n\n";
      }
      elsif (-e "$curdir/#DONE#")
      {
        print "Process ended on $RunningProcesses{$curdir} at ",scalar localtime(time()),
        " (", (scalar @cpulist)," processors available).\n\n";
      }
      else
      {
        # unexpected case.  Some sort of strangeness is going on.
        print "Process exited with UNKNOWN STRANGENESS on $RunningProcesses{$curdir} at ",
        scalar localtime(time()), " (", (scalar @cpulist)," processors available).\n\n";
      }
    }
    else
    {
      # add to list of still running jobs
      push @remainingjobs, $_;
    }

  }

  # update running job list
  @runningjobs = @remainingjobs;
}

#--------------------------------------------------------------------------
sub start_waiting_jobs
{
  #  uses globals:
  #     cpulist
  #     RootBFBDir
  #     runningjobs
  #     RunningProcesses
  #
  if (scalar @cpulist > 0)  ####### if a processor is ready get it started
  {
    ###### one or more processors are available, so look for a path to compute
    my @pths = ();
    find_names($RootBFBDir,"#$PROCESSFLAG#",\@pths);
    foreach (@pths)
    {
      my $found = $_;
      if (scalar @cpulist > 0)
      {
        my $cpu = shift @cpulist;
        
        my $newdir = $found;
        $newdir =~ s/\/#$PROCESSFLAG#//;
        unlink($found);
        
        my @tmp = split('/',$newdir);
        my $flnm = pop @tmp;
        
        # update bfb and in files
        system("gunzip -f $newdir/$flnm.in.gz >& /dev/null");
        
        # run it
        my $retval = system("ssh $cpu \"(cd $newdir; touch \\#PROCESSING#; $ProgExec < $flnm.in >& $flnm.out; if [ -e abort.dat ]; then /bin/mv -f \\#PROCESSING# \\#ABORTED#; elif [ \\`grep --count 'QC simulation terminated' qc.log\\` == 1 ]; then /bin/mv -f \\#PROCESSING# \\#DONE#; else /bin/mv -f \\#PROCESSING# \\#ERROR#; fi) < /dev/null >& /dev/null &\"");
        
        # update runningjobs list
        push @runningjobs, $newdir;
        # update hash of running processes hash
        $RunningProcesses{$newdir} = $cpu;
        # print message to screen
        print "Process started on $cpu to compute $flnm at ",scalar localtime(time()),
        " (",(scalar @cpulist)," processors available).\n\n";
      }
    }
  }
}

#--------------------------------------------------------------------------
sub find_names
{
  # Arg list
  #  - root directory name
  #  - string to look for
  #  - list reference for found items to be pushed onto
  
  my $dirname = shift;
  my $searchstr = shift;
  my $listref = shift;
  my @lst;
  
  $wantedref = sub {@lst = split('/',$_); $_ = pop @lst;
    if ($_ =~ $searchstr)
    {
      push @{$listref}, $File::Find::name;
    }
  };
  
  find({wanted => $wantedref, no_chdir => 1}, $dirname);
}

# finds and returns first file matching
sub find_first
{
  # Arg list
  #  - root directory name
  #  - string to look for

  my $dirname = shift;
  my $searchstr = shift;
  my @lst;

  find_names($dirname,$searchstr,\@lst);

  return $lst[0];
}


#################################################################################

sub ceiling
{
  my $val = shift;
  my $bot = int($val);

  if ($val > $bot)
  {
    return ($bot+1);
  }
  else
  {
    return $bot;
  }
}
