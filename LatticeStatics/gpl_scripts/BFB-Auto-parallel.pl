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
        "BFB-auto <wall-clock time in seconds for when to START shutting down> <rel-path-name-of-executable> <rel-path-to-bfbrootdir>\n";
}

$walltimeout = shift;
$ProgExec = getcwd() . "/" . shift;
$RootBFBDir = getcwd() . "/" . shift;

STDIN->autoflush(1);
print "\n\nStarting Parallel QCBFB\n\n";

# global variables
$InputFileName = "columnNi";
$sym_info = "symmetry_bases.pm";
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


############ clean up aborted directories and files ###############
print "Processing ABORTED directories.\n";
$i=0;


@abtd = ();
find_names($RootBFBDir,"#ABORTED#",\@abtd);
foreach (@abtd)
{
  $i++;
  $workingdir = $_;
  $workingdir =~ s/\/#ABORTED#//;

  @dummy = split('/',$workingdir); 
  $shortdir = pop @dummy;

  print "     Working on: ", $shortdir,"\n";
  
  if ($workingdir eq $RootBFBDir)
  {
    system("gunzip -f $workingdir/$InputFileName.bfb.gz $workingdir/$InputFileName.in.gz >& /dev/null");
  }
  else
  {
    system("gunzip -f $workingdir/$shortdir.bfb.gz $workingdir/$shortdir.in.gz $workingdir/$shortdir.res.gz >& /dev/null");
  }
  unlink glob("$workingdir/*.plt.gz");
  unlink "$workingdir/abort.dat";
  unlink glob("$workingdir/*.out.gz");
  unlink glob("$workingdir/*.B[0-9]*.gz");
  unlink glob("$workingdir/*.T[0-9]*.gz");
  unlink glob("$workingdir/*.E[0-9]*.gz");
  unlink glob("$workingdir/*.bpp.gz");
  unlink "$workingdir/qc.log.gz";
  unlink "$workingdir/qc.cmd.gz";

  move("$workingdir/#ABORTED#","$workingdir/#WAITING#");
}
if ($i == 0)
{
  print "     No ABORTED directories found.\n";
}
print "Done processing ABORTED direcotries.\n\n";


################## find number of point on path requested #########################
if (-e "$RootBFBDir/$InputFileName.in")
{
   open(INFL,"$RootBFBDir/$InputFileName.in") or die "can't open input file to find num. pts.\n";
}
else
{
   $zipped = 1;
   system("gunzip -f $RootBFBDir/$InputFileName.in.gz >& /dev/null");
   open(INFL,"$RootBFBDir/$InputFileName.in") or die "can't open input file to find num. pts.\n";
}
$_ = <INFL>;
while (defined($_) && ($_ !~ /^bfb/))
{
   $_ = <INFL>;
}
while (defined($_) && ($_ !~ /^loop/))
{
   $_ = <INFL>;
}
close(INFL);
if ($zipped == 1)
{
  system("gzip -f $RootBFBDir/$InputFileName.in >& /dev/null");
}
@t = split ',', $_;
$NumPts = pop @t;
chomp $NumPts;

print "Using $NumPts points on each BFB curve.\n\n";


############## spawn process to look for sentinel files ##################
print "Spawning sentinel process... ";

open(SENTINELREPORT,">$RootBFBDir/SentinelReport.txt");
SENTINELREPORT->autoflush(1);
defined($SentinelPID = fork()) or die "can't fork sentinel.\n";
if ($SentinelPID == 0)
{
  # sentinel's job
  while (-e $maintimerfile)
  {
    find_sentinels($RootBFBDir,$NumPts);

    print SENTINELREPORT ".";
    sleep($sleeptime);
    print SENTINELREPORT ".";
  }
  
  print SENTINELREPORT "Exiting ", scalar localtime(time()), "   .\n";
  exit;
}
print "Done. Sentinel pid = $SentinelPID.\n\n";


####################### start running #####################
# check root
if (! ( (-e "$RootBFBDir/#DONE#") || (-e "$RootBFBDir/#ERROR#") ))
{
  open(ROT,">$RootBFBDir/#WAITING#");
  close(ROT);
  if (-e "$RootBFBDir/$InputFileName.in")
  {
    system("gzip -f $RootBFBDir/$InputFileName.in >& /dev/null");
  }
  if (-e "$RootBFBDir/$InputFileName.bfb")
  {
    system("gzip -f $RootBFBDir/$InputFileName.bfb >& /dev/null");
  }
}

while((-e $maintimerfile) &&
      ( ((scalar @cpulist) < $TotalNumProcessors) ||
        (defined(find_first($RootBFBDir,"#WAITING#"))) ||
        (defined(find_first($RootBFBDir,"sentinel")))
      )
     )
{
  #### find any processors that have finished
  @pths = ();
  find_names($RootBFBDir,"#RUNNING#",\@pths);
  foreach (@pths)
  {
    $curdir = $_;
    $curdir =~ s/\/#RUNNING#//;
    @tmp = split('/',$curdir);
    $flnm = pop @tmp;

    if (! -e "$curdir/#PROCESSING#")
    {
      #### process has finished
      # clean up after the run returns.
      
      # use -q to stop warnings
      if ($newdir ne $RootBFBDir)
      {
        # if not the root path
        system("gzip -f $curdir/$flnm.bfb $curdir/$flnm.in $curdir/$flnm.res $curdir/$flnm.out $curdir/*.plt" .
           " $curdir/${flnm}curdir*.res $curdir/*.T[0-9]*.bfb $curdir/*.T[0-9]*.res $curdir/*.E[0-9]*.bfb $curdir/*.E[0-9]*.res $curdir/$flnm.bpp $curdir/qc.* >& /dev/null");
      }
      else
      {
        # if the root path
        # set the input file root
        $flnm = $InputFileName;
        system("gzip -f $curdir/$flnm.bfb $curdir/$flnm.in $curdir/$flnm.res $curdir/$flnm.out $curdir/*.plt" .
           " $curdir/${flnm}curdir*.res $curdir/*.T[0-9]*.bfb $curdir/*.T[0-9]*.res $curdir/*.E[0-9]*.bfb $curdir/*.E[0-9]*.res $curdir/$flnm.bpp $curdir/qc.* >& /dev/null");
      }

      if (-e "$newdir/abort.dat")
      {
        move($curdir . "/#RUNNING#", $curdir . "/#ABORTED#");
      }
      else
      {
        move($curdir . "/#RUNNING#", $curdir . "/#DONE#");
      }
      push @cpulist, $RunningProcesses{$curdir};
      print "Process ended on $RunningProcesses{$curdir} at ",scalar localtime(time()),
            " (", (scalar @cpulist)," processors available).\n\n";
      delete $RunningProcesses{$curdir};
    }
    elsif ( abs(( (stat("$curdir/#PROCESSING#"))[9] - time() )/60.0) > 5.0 ) # if older than 5 minutes
    {
      #### process has exited
      
      # use -q to stop warnings
      if ($newdir ne $RootBFBDir)
      {
        # if not the root path
        system("gzip -f $curdir/$flnm.bfb $curdir/$flnm.in $curdir/$flnm.res $curdir/$flnm.out $curdir/*.plt" .
           " $curdir/${flnm}curdir*.res $curdir/*.T[0-9]*.bfb $curdir/*.T[0-9]*.res $curdir/*.E[0-9]*.bfb $curdir/*.E[0-9]*.res $curdir/$flnm.bpp $curdir/qc.* >& /dev/null");
      }
      else
      {
        # if the root path
        # set the input file root
        $flnm = $InputFileName;
        system("gzip -f $curdir/$flnm.bfb $curdir/$flnm.in $curdir/$flnm.res $curdir/$flnm.out $curdir/*.plt" .
           " $curdir/${flnm}curdir*.res $curdir/*.T[0-9]*.bfb $curdir/*.T[0-9]*.res $curdir/*.E[0-9]*.bfb $curdir/*.E[0-9]*.res $curdir/$flnm.bpp $curdir/qc.* >& /dev/null");
      }

      move($curdir . "/#RUNNING#", $curdir . "/#ERROR#");
      unlink("$curdir/#PROCESSING#");
      push @cpulist, $RunningProcesses{$curdir};
      print "Process exited (ERROR) on $RunningProcesses{$curdir} at ",scalar localtime(time()),
            " (", (scalar @cpulist)," processors available).\n\n";
      delete $RunningProcesses{$curdir};
    }
  }

  if (scalar @cpulist > 0)  ####### if a processor is ready get it started
  {
    ###### one or more processors are available, so look for a path to compute
    @pths = ();
    find_names($RootBFBDir,"#WAITING#",\@pths);
    foreach (@pths)
    {
      $found = $_;
      if (scalar @cpulist > 0)
      {
        $cpu = shift @cpulist;
        
        $newname = $found;
        $newname =~ s/WAITING/RUNNING/;
        move($found,$newname);
        
        $newdir = $newname;
        $newdir =~ s/\/#RUNNING#//;
        @tmp = split('/',$newdir);
        $flnm = pop @tmp;
        
        if ($newdir ne $RootBFBDir)
        {
          # if not the root path
          # update bfb and in files
          system("gunzip -f $newdir/$flnm.bfb.gz $newdir/$flnm.in.gz $newdir/$flnm.res.gz >& /dev/null");
        }
        else
        {
          # if the root path
          # set the input file root
          $flnm = $InputFileName;
          system("gunzip -f $newdir/$flnm.in.gz $newdir/$flnm.bfb.gz >& /dev/null");
        }
        
        # run it
        $retval = system("ssh $cpu \"(cd $newdir; touch \\#PROCESSING#; $ProgExec < $flnm.in >& $flnm.out; if [ \\`grep --count 'QC simulation terminated' qc.log\\` == 1 ]; then /bin/rm -f \\#PROCESSING#; fi) < /dev/null >& /dev/null &\"");
        
        # update hash of running processes
        print "Process started on $cpu to compute $flnm at ",scalar localtime(time()),
        " (",(scalar @cpulist)," processors available).\n\n";
        $RunningProcesses{$newdir} = $cpu;
      }
    }
  }
  
  if (scalar @cpulist > 0)
  {
    # there are no paths ready to be computed then sleep
    sleep($sleeptime);
  }
}
# finished main loop

# set abort.dat for any running processes
@stillrunning = ();
find_names($RootBFBDir,"#RUNNING#",\@stillrunning);
foreach(@stillrunning)
{
  $abortdir = $_;
  $abortdir =~ s/\/#RUNNING#//;

  $flnm = $abortdir;
  $flnm =~ s/.*\/([^\/]*)$/$1/;
  print "Aborting $flnm.\n";
  open(ABORTFL,">$abortdir/abort.dat");
  print ABORTFL "T\n";
  close(ABORTFL);
}
print "\n";

# check it time remains
if (-e "$maintimerfile")
{
  print "Finished and Timer file still exists.  Removing.\n";
  kill 9, $TimerPID;
  unlink "$maintimerfile";
}

# wait for all child processes to exit
print "Waiting for child processes to end...\n";
#### find any processors that have finished
@pths = ();
find_names($RootBFBDir,"#RUNNING#",\@pths);
while ( (scalar @pths) > 0)
{
  foreach (@pths)
  {
    $curdir = $_;
    $curdir =~ s/\/#RUNNING#//;
    @tmp = split('/',$curdir);
    $flnm = pop @tmp;

    if (! -e "$curdir/#PROCESSING#")
    {
      #### process has finished
      # clean up after the run returns.
      
      # use -q to stop warnings
      if ($newdir ne $RootBFBDir)
      {
        # if not the root path
        system("gzip -f $curdir/$flnm.bfb $curdir/$flnm.in $curdir/$flnm.res $curdir/$flnm.out $curdir/*.plt" .
           " $curdir/${flnm}curdir*.res $curdir/*.T[0-9]*.bfb $curdir/*.T[0-9]*.res $curdir/*.E[0-9]*.bfb $curdir/*.E[0-9]*.res $curdir/$flnm.bpp $curdir/qc.* >& /dev/null");
      }
      else
      {
        # if the root path
        # set the input file root
        $flnm = $InputFileName;
        system("gzip -f $curdir/$flnm.bfb $curdir/$flnm.in $curdir/$flnm.res $curdir/$flnm.out $curdir/*.plt" .
           " $curdir/${flnm}curdir*.res $curdir/*.T[0-9]*.bfb $curdir/*.T[0-9]*.res $curdir/*.E[0-9]*.bfb $curdir/*.E[0-9]*.res $curdir/$flnm.bpp $curdir/qc.* >& /dev/null");
      }
      
      if (-e "$newdir/abort.dat")
      {
        move($curdir . "/#RUNNING#", $curdir . "/#ABORTED#");
      }
      else
      {
        move($curdir . "/#RUNNING#", $curdir . "/#DONE#");
      }
      push @cpulist, $RunningProcesses{$curdir};
      print "     Process ended on $RunningProcesses{$curdir} ",
      "(", (scalar @cpulist)," processors available).\n\n";
      delete $RunningProcesses{$curdir};
    }
    elsif ( abs(( (stat("$curdir/#PROCESSING#"))[9] - time() )/60.0) > 5.0 ) # if older than 5 minutes
    {
      #### process has exited

      # use -q to stop warnings
      if ($newdir ne $RootBFBDir)
      {
        # if not the root path
        system("gzip -f $curdir/$flnm.bfb $curdir/$flnm.in $curdir/$flnm.res $curdir/$flnm.out $curdir/*.plt" .
           " $curdir/${flnm}curdir*.res $curdir/*.T[0-9]*.bfb $curdir/*.T[0-9]*.res $curdir/*.E[0-9]*.bfb $curdir/*.E[0-9]*.res $curdir/$flnm.bpp $curdir/qc.* >& /dev/null");
      }
      else
      {
        # if the root path
        # set the input file root
        $flnm = $InputFileName;
        system("gzip -f $curdir/$flnm.bfb $curdir/$flnm.in $curdir/$flnm.res $curdir/$flnm.out $curdir/*.plt" .
           " $curdir/${flnm}curdir*.res $curdir/*.T[0-9]*.bfb $curdir/*.T[0-9]*.res $curdir/*.E[0-9]*.bfb $curdir/*.E[0-9]*.res $curdir/$flnm.bpp $curdir/qc.* >& /dev/null");
      }

      move($curdir . "/#RUNNING#", $curdir . "/#ERROR#");
      unlink("$curdir/#PROCESSING#");
      push @cpulist, $RunningProcesses{$curdir};
      print "     Process exited (ERROR) on $RunningProcesses{$curdir} ",
      "(", (scalar @cpulist)," processors available).\n\n";
      delete $RunningProcesses{$curdir};
    }
  }

  sleep($sleeptime);
  @pths = ();
  find_names($RootBFBDir,"#RUNNING#",\@pths);
}

while (-1 != ($ans =wait))
{
  print "     Process $ans ended.\n";
}
print "Done.\n\n";

# run final sentinel check
find_sentinels($RootBFBDir,$NumPts);

# clean up
close(SENTINELREPORT);

print "Ending Parallel QCBFB.\n\n";
exit;

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

# find any waiting sentinel files
#
# argument list
#
# - root directory to search
#
sub find_sentinels
{
  my $rootdir = shift;
  my $NumPts = shift;
  my $dir;
  my $found;
  my @sentls;

  @sentls = ();
  find_names($rootdir,"sentinel",\@sentls);
  foreach (@sentls)
  {
    $found = $_;
    $dir = $found;
    $dir =~ s/(.*)\/[^\/]*/$1/;
    $found =~ s/.*\/([^\/]*)\.sentinel/$1/;
   
    print SENTINELREPORT "\nProcessing $found...";
    if (!-d "$dir/$found")
    { 
      find_sym_and_update_bfb("$dir","$found.bfb");
      update_in_file("$dir/$found.in",$NumPts);

      mkdir "$dir/$found";
      foreach (glob("$dir/$found.*"))
      {
        if (!/sentinel/)
        {
          system("gzip -f $_ >& /dev/null");
          move("$_.gz", "$dir/$found/");
        }
      }
      symlink "../$sym_info", "$dir/$found/$sym_info";
      symlink "../$geo","$dir/$found/$geo";
      symlink "../$pots","$dir/$found/$pots";

      if ($found =~ /\.B....-...\.B....-...\.B....-...\.B....-.../)
      {
        system("gzip -f $dir/$found/$found.bfb $dir/$found/$found.in $dir/$found/$found.res >& /dev/null");
        open(WAT,">$dir/$found/#SKIPPED#");
        close(WAT);
      }
      else
      {
        open(WAT,">$dir/$found/#WAITING#");
        close(WAT);
      }
      unlink "$dir/$found.sentinel";
    }
    else
    {
      foreach (glob("$dir/$found.*"))
      {
        unlink $_;
      }
    }
    print SENTINELREPORT "   Done at ",scalar localtime(time()), ". ";
  }
}

sub find_sym_and_update_bfb
{
  my $curdir = shift;
  my $flnm = shift;

  do "$curdir/symmetry_bases.pm";
  
  my $olddir = cwd();
  chdir $curdir;
  do "$curdir/$flnm";
  chdir $olddir;

  my ($T,$Tol,$a,$aa,$aaa,$b,$bb,$bbb,$c,$cc,$ccc,$d,$dd,$ddd,$e,$ee,$eee,$f,$ff,$fff,$g,$gg,
      $ggg,$IsTRx,$IsTRy,$IsTRz,$IsTQx,$IsTQy,$IsTQz,$IsTJ,$SymGrp,$fl,$foundsymmat,$dummy,$type);
      
  @T = @{$StartType{Tangent}};

  $Tol = $SolutionMethod{NewtonPCSolution}{ConvergeCriteria};
  pop @T;
  
  $a = matrix_product(\@TRx,\@T);
  $aa = difference($a,\@T);
  $aaa = dot_product($aa,$aa);
  $IsTRx = ($aaa <= $Tol);
  #print "norm(TRx*Tangent - Tangent) = ",$aaa,"\n";
  
  $b = matrix_product(\@TRy,\@T);
  $bb = difference($b,\@T);
  $bbb = dot_product($bb,$bb);
  $IsTRy = ($bbb <= $Tol);
  #print "norm(TRy*Tangent - Tangent) = ",$bbb,"\n";

  $c = matrix_product(\@TRz,\@T);
  $cc = difference($c,\@T);
  $ccc = dot_product($cc,$cc);
  $IsTRz = ($ccc <= $Tol);
  #print "norm(TRz*Tangent - Tangent) = ",$ccc,"\n";

  $d = matrix_product(\@TQx,\@T);
  $dd = difference($d,\@T);
  $ddd = dot_product($dd,$dd);
  $IsTQx = ($ddd <= $Tol);
  #print "norm(TQx*Tangent - Tangent) = ",$ddd,"\n";

  $e = matrix_product(\@TQy,\@T);
  $ee = difference($e,\@T);
  $eee = dot_product($ee,$ee);
  $IsTQy = ($eee <= $Tol);
  #print "norm(TQy*Tangent - Tangent) = ",$eee,"\n";

  $f = matrix_product(\@TQz,\@T);
  $ff = difference($f,\@T);
  $fff = dot_product($ff,$ff);
  $IsTQz = ($fff <= $Tol);
  #print "norm(TQz*Tangent - Tangent) = ",$fff,"\n";

  $g = matrix_product(\@TJ,\@T);
  $gg = difference($g,\@T);
  $ggg = dot_product($gg,$gg);
  $IsTJ = ($ggg <= $Tol);
  #print "norm(TJ*Tangent - Tangent) = ",$ggg,"\n";

  if ($IsTRx && $IsTRy && $IsTRz && $IsTQx && $IsTQy && $IsTQz && $IsTJ)
  {
    $SymGrp = "G0";
    @SymChkList = ();
  }
  elsif ($IsTRx && $IsTRy && $IsTQz)
  {
    $SymGrp = "G4";
    @SymChkList = ("C0_4");
  }
  elsif ($IsTRx && $IsTRz && $IsTQy)
  {
    $SymGrp = "G5";
    @SymChkList = ("C0_5");
  }
  elsif ($IsTRy && $IsTRz && $IsTQx)
  {
    $SymGrp = "G6";
    @SymChkList = ("C0_6");
  }
  elsif ($IsTRz && $IsTQz && $IsTJ)
  {
    $SymGrp = "G7";
    @SymChkList = ("C0_7");
  }
  elsif ($IsTRx && $IsTQx && $IsTJ)
  {
    $SymGrp = "G12";
    @SymChkList = ("C0_12");
  }
  elsif ($IsTRy && $IsTQy && $IsTJ)
  {
    $SymGrp = "G13";
    @SymChkList = ("C0_13");
  }
  elsif ($IsTQx && $IsTQy && $IsTQz)
  {
    $SymGrp = "G14";
    @SymChkList = ("C0_14");
  }
  elsif ($IsTRx)
  {
    $SymGrp = "G1";
    @SymChkList = ("C4_1", "C5_1", "C12_1");
  }
  elsif ($IsTRy)
  {
    $SymGrp = "G2";
    @SymChkList = ("C4_2", "C6_2", "C13_2");
  }
  elsif ($IsTRz)
  {
    $SymGrp = "G3";
    @SymChkList = ("C5_3", "C6_3", "C7_3");
  }
  elsif ($IsTQx)
  {
    $SymGrp = "G8";
    @SymChkList = ("C6_8", "C12_8", "C14_8");
  }
  elsif ($IsTQy)
  {
    $SymGrp = "G9";
    @SymChkList = ("C5_9", "C13_9", "C14_9");
  }
  elsif ($IsTQz)
  {
    $SymGrp = "G10";
    @SymChkList = ("C4_10", "C7_10", "C14_10");
  }
  elsif ($IsTJ)
  {
    $SymGrp = "G11";
    @SymChkList = ("C7_11", "C12_11", "C13_11");
  }
  else
  {
    $SymGrp = "Id"; # G15
    @SymChkList = ("C1_15", "C2_15", "C3_15", "C8_15", "C9_15", "C10_15", "C11_15");
  }

  #print "has symmetry group $SymGrp\n";
  
  $fl='';
  $tangent='';
  $bifpt='';
 
 
  $fl .= "use symmetry_bases;\n";  # put this at the top to make sure it is loaded before it is needed
  open(ORIGFL,"$curdir/$flnm");
  while (<ORIGFL>)
  {
    # assume tangent and bifpt are one line.  keep only the last one found.
    if (/StartType{Tangent}/)
    {
      $tangent = $_;
    }
    elsif (/StartType{BifurcationPoint}/)
    {
      $bifpt = $_;
    }
    elsif (/Restriction{Type}/)
    {
      if ($SymGrp eq "Id")
      {
        $fl .= "\$Restriction{Type} = NoRestriction;\n";
      }
      else
      {
        $fl .= "\$Restriction{Type} = RestrictToTranslatedSubSpace;\n";
      }

      $fl .= "\$Restriction{SymmetryCheckProjectionMatrices} = [";
      $fl .= "[@" . (shift @SymChkList) . "]";
      foreach $mat (@SymChkList)
      {
        $fl .= ",[@" . $mat . "]";
      } 
      $fl .= "];\n";
    }
    elsif (/Restriction{SymmetryCheckProjectionMatrices}/)
    {
      # Remove all symmetry matrices
      while (defined($_) && ($_ !~ /.*];$/))
      {
        $_=<ORIGFL>;
      }
    }
    elsif (/{RestrictToTranslatedSubSpace}{ProjectionMatrix}/)
    {
      if ($SymGrp ne "Id")
      {
        $fl .= "\$Restriction{RestrictToTranslatedSubSpace}{ProjectionMatrix} = [@" 
            . $SymGrp . "];\n";
      }

      while (defined($_) && ($_ !~ /.*];$/))
      {
        $_=<ORIGFL>;
      }
    }
    else
    {
      $fl .= $_;
    }
  }
  close(ORIGFL);

  # put in the tangent and the bif pt
  $fl .= $tangent;
  $fl .= $bifpt;

  unlink("$curdir/$flnm");
  open(NEWFL,">$curdir/$flnm");
  printf NEWFL "$fl";
  close(NEWFL);
}

sub update_in_file
{
  # turn on inplace file editing
  my $NumPts = pop @_;
  open INFILE, "$_[0]";
  open OUTFILE, ">$_[0]~";
  while (<INFILE>)
  {
    if (/^loop/)
    {
       print OUTFILE "loop,,$NumPts\n";
    }
    else
    {
      print OUTFILE $_;
    }
  }
  close INFILE;
  close OUTFILE;
  move("$_[0]~", "$_[0]");
}

#################################################################################

# arguments T,v -- S = T*v
sub matrix_product
{
  my($T) = shift;
  my(@TT)=@{$T};
  my($v) = shift;
  my(@vv) = @{$v};
  my($i,$j,$len);
  my(@S);

  $len = scalar @vv;
  for ($i=0;$i<$len;$i++)
  {
    $S[$i] = 0.0;
    for ($j=0;$j<$len;$j++)
    {
      $S[$i] += $TT[$i][$j]*$vv[$j];
    }
  }

  return \@S;
}

sub dot_product
{
   my($v) = shift;
   my(@vv) = @{$v};
   my($w) = shift;
   my(@ww) = @{$w};
   my($val,$i,$len);

   $val = 0.0;
   $len = scalar @vv;
   for ($i=0;$i<$len;$i++)
   {
     $val += $vv[$i]*$ww[$i];
   }

   return $val;
}

sub difference
{
  my($v) = shift;
   my(@vv) = @{$v};
   my($w) = shift;
   my(@ww) = @{$w};
   my(@val,$i,$len);

   $len = scalar @vv;
   for ($i=0;$i<$len;$i++)
   {
      $val[$i] = $vv[$i] - $ww[$i];
   }

   return \@val;
}

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
