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

if ((scalar @ARGV) == 0)
{
  die "Usage\n" .
        "BFB-auto <wall-clock time in seconds for when to STAET shutting down> <rel-path-name-of-executable> <rel-path-to-bfbrootdir>\n";
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
    $cpulist[$i] = $_;
    print "Processor $i has name: ",$cpulist[$i],"\n";
    $i++;
    $TotalNumProcessors++;
  }
}
close(LST);
print "Total number of processors is $TotalNumProcessors.\n\n";
############ create runtime sentinel and timer ######################

$maintimerfile = "$RootBFBDir/#MASTER-TIMER#";
print "Seting wall-clock timeout to $walltimeout seconds (",$walltimeout/60," minutes).... ";
defined($TimerPID = fork()) or die "can't fork timer.\n";
if ($TimerPID == 0)
{
  # Timer's job
  $remainingtime = $walltimeout;
  if ($decrement > $remainingtime)
  {
    $decrement = $remainingtime;
  }
  open(MAINTIMEFILE,">$maintimerfile");
  MAINTIMEFILE->autoflush(1);

  print MAINTIMEFILE "$remainingtime seconds (",$remainingtime/60," minutes) remaining at ", scalar localtime(time()), ".\n";
  while($remainingtime > 0)
  {
    sleep($decrement);
    $remainingtime -= $decrement;
    print MAINTIMEFILE "$remainingtime seconds (",$remainingtime/60," minutes) remaining at ", scalar localtime(time()), ".\n";
  }

  unlink("$maintimerfile");

  exit;
}
print "Done. Timer pid = $TimerPID.\n\n";


############ clean up aborted directories and files ###############
print "Processing ABORTED directories.\n";
$i=0;
while (defined($_ = find_first($RootBFBDir,"#ABORTED#")))
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
  unlink glob("$workingdir/*BP*.gz");
  unlink glob("$workingdir/*TP*.gz");
  unlink glob("$workingdir/*.bpp.gz");
  unlink "$workingdir/qc.log.gz";

  move("$workingdir/#ABORTED#","$workingdir/#WAITING#");
}
if ($i == 0)
{
  print "     No ABORTED directories found.\n";
}
print "Done processing ABORTED direcotries.\n\n";


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
    find_sentinels($RootBFBDir);

    sleep($sleeptime);
  }
  
  exit;
}
print "Done. Sentinel pid = $SentinelPID.\n\n";


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
while ($_ !~ /^bfb/)
{
   $_ = <INFL>;
}
while ($_ !~ /^loop/)
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


####################### start running #####################
# check root
if (! -e "$RootBFBDir/#DONE#")
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
  #### if no processors available wait for next to finish
  if ((scalar @cpulist) == 0)
  {
    print "Waiting for next process to end...\n";
    $pidend = wait;
    if ($pidend == $SentinelPID)
    {
      print "Sentinel process ended.\n\n";
    }
    elsif ($pidend == $TimerPID)
    {
      print "Timer process ended.\n\n";
    }
    else
    {
      print "Process $pidend ended on $RunningProcesses{$pidend} (1 processor available).\n\n";
      push @cpulist, $RunningProcesses{$pidend};
      delete $RunningProcesses{$pidend};
    }
  }
  else ###### else look for other processes that have ended
  {
    foreach (keys %RunningProcesses)
    {
      $ans = waitpid($_, &WNOHANG);
      if ($ans != 0)
      {
        push @cpulist, $RunningProcesses{$ans};
        print "Process $ans ended on $RunningProcesses{$ans} ",
              "(",(scalar @cpulist)," processors available).\n\n";
        delete $RunningProcesses{$ans};
      }
    }
  }
  
  if ((scalar @cpulist) > 0)
  {
    ###### a processor is available, so look for a path to compute
    $found = find_first($RootBFBDir,"#WAITING#");
    if (defined($found))  ####### if a path is ready get it started
    {
      $cpu = shift @cpulist;

      # need to move to #RUNNING# before fork to ensure that it is not started multiple times
      $newname = $found;
      $newname =~ s/WAITING/RUNNING/;
      move($found,$newname);
      if (!defined($pid = fork()))
      {
        move($newname,$found);
        die "can't fork child.\n";
      }
      if ($pid == 0)
      { ##### this is the child
        $newdir = $newname;
        $newdir =~ s/\/#RUNNING#//;
        chdir($newdir);
        @tmp = split('/',$newdir);
        $flnm = pop @tmp;

        
        if ($newdir ne $RootBFBDir)
        {
          # if not the root path
          # update bfb and in files
          system("gunzip -f $flnm.bfb.gz $flnm.in.gz $flnm.res.gz >& /dev/null");
          find_sym_and_update_bfb("$flnm.bfb");
          update_in_file("$flnm.in",$NumPts);
        }
        else
        {
          # if the root path
          # set the input file root
          $flnm = $InputFileName;
          system("gunzip -f $flnm.in.gz $flnm.bfb.gz >& /dev/null");
        }

        # run it
        $retval = system("ssh $cpu \"(cd $newdir && $ProgExec < $flnm.in >& $flnm.out;  sync)\"");

        # clean up after the run.
      
        # use -q to stop warnings
        system("gzip -f $newdir/$flnm.bfb $newdir/$flnm.in $newdir/$flnm.res $newdir/$flnm.out $newdir/*.plt" .
               " $newdir/${flnm}_*.res $newdir/*TP*.bfb $newdir/*TP*.res $newdir/$flnm.bpp $newdir/qc.* >& /dev/null");
        if ($retval != 0) # error
        {
          move($newdir . "/#RUNNING#", $newdir . "/#ERROR#");
        }
        elsif (-e "$newdir/abort.dat")
        {
          move($newdir . "/#RUNNING#", $newdir . "/#ABORTED#");
        }
        else
        {
          move($newdir . "/#RUNNING#", $newdir . "/#DONE#");
        }
        exit;
        ##### end of child
      }
    
      # parent continues and updates hash of running processes
      $newdir = $found;
      $newdir =~ s/\/#WAITING#//;
      @tmp = split('/',$newdir);
      $flnm = pop @tmp;
      print "Process $pid started on $cpu to compute $flnm at ",scalar localtime(time()), ".\n\n";
      $RunningProcesses{$pid} = $cpu;
    }
    else
    {
      # if no paths are ready to be computed then sleep
      sleep($sleeptime);
    }
  }
}
# finished main loop

# set abort.dat for any running processes
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
while (-1 != ($ans =wait))
{
  if (($ans != $TimerPID) && ($ans != $SentinelPID))
  {
    push @cpulist, $RunningProcesses{$ans};
    print "    Process $ans ended on $RunningProcesses{$ans} ",
          "(",(scalar @cpulist)," processors available).\n";
    delete $RunningProcesses{$ans};
  }
  else
  {
    print "    Process $ans ended.\n";
  }
}
print "Done.\n\n";

# run final sentinel check
find_sentinels($RootBFBDir);

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
  my $dir;
  my $found;

  while(defined($found = find_first($rootdir,"sentinel")))
  {
    $dir = $found;
    $dir =~ s/(.*)\/[^\/]*/$1/;
    $found =~ s/.*\/([^\/]*)\.sentinel/$1/;
   
    print SENTINELREPORT "Processing $found...";
    if (!-d "$dir/$found")
    { 
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
      if ($found =~ /\.BP\....\.BP\....\.BP\..../)
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
    print SENTINELREPORT "   Done at ",scalar localtime(time()), ".\n";
  }
}

sub find_sym_and_update_bfb
{
  my $flnm = shift;
  do "symmetry_bases.pm";
  do "$flnm";

  my ($T,$Tol,$a,$aa,$aaa,$b,$bb,$bbb,$c,$cc,$ccc,$d,$dd,$ddd,$e,$ee,$eee,$f,$ff,$fff,$g,$gg,
      $ggg,$IsTRx,$IsTRy,$IsTRz,$IsTQx,$IsTQy,$IsTQz,$IsTJ,$SymGrp,$fl);
      

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
  }
  elsif ($IsTRx && $IsTRy && $IsTQz)
  {
    $SymGrp = "G4";
  }
  elsif ($IsTRx && $IsTRz && $IsTQy)
  {
    $SymGrp = "G5";
  }
  elsif ($IsTRy && $IsTRz && $IsTQx)
  {
    $SymGrp = "G6";
  }
  elsif ($IsTRz && $IsTQz && $IsTJ)
  {
    $SymGrp = "G7";
  }
  elsif ($IsTRx && $IsTQx && $IsTJ)
  {
    $SymGrp = "G12";
  }
  elsif ($IsTRy && $IsTQy && $IsTJ)
  {
    $SymGrp = "G13";
  }
  elsif ($IsTQx && $IsTQy && $IsTQz)
  {
    $SymGrp = "G14";
  }
  elsif ($IsTRx)
  {
    $SymGrp = "G1";
  }
  elsif ($IsTRy)
  {
    $SymGrp = "G2";
  }
  elsif ($IsTRz)
  {
    $SymGrp = "G3";
  }
  elsif ($IsTQx)
  {
    $SymGrp = "G8";
  }
  elsif ($IsTQy)
  {
    $SymGrp = "G9";
  }
  elsif ($IsTQz)
  {
    $SymGrp = "G10";
  }
  elsif ($IsTJ)
  {
    $SymGrp = "G11";
  }
  else
  {
    die "Can't determine subgroup...\n";
  }

  #print "has symmetry group $SymGrp\n";
  
  $fl='';
  
  open(ORIGFL,$flnm);
  while (<ORIGFL>)
  {
    if (/Restriction{Type}/)
    {
      $fl .= $_;
      $_ = <ORIGFL>;
      if (!/symmetry_bases/)
      {
      $fl .= "use symmetry_bases;\n";
      }
    }
    
    if (/{RestrictToTranslatedSubSpace}{ProjectionMatrix}/)
    {
      $fl .= "\$Restriction{RestrictToTranslatedSubSpace}{ProjectionMatrix} = [@" 
          . $SymGrp . "];\n";
      while ($_ !~ /.*];$/)
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
  
  open(NEWFL,">$flnm");
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
