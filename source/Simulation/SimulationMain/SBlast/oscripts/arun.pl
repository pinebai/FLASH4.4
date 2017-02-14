#!/usr/bin/perl -w

# Run the Sedov problem with FLASH over the
# parameter space: source energy, source density

my $title    = "Sedov";
my $term     = "png";
my $dotterm  = ".png";
my $lw       = 2;

my $dekfile  = "flash1d_uni.par";                  # Prototype flash.par file
my $parfile  = "this.par";                         # flash.par file used for run

my $idl      = "/usr/local/rsi/idl_6.1/bin/idl";   # Path to idl
my $sedovsolver = "dumps/sedov";                   # Path to sedov analytical solver

my $max_dumps = 100;                               # Maximum number of chekpoint files

my @colors = (4,2,3,1);                            # gnuplot line styles 

my @EIn = (1);                                     # Source energy
my @rhoIn = (1);                                   # Source density

my %sim_params;

my @r=(); 
my @rho=();
my @p=(); 
my @u=();

my @rho_sol=();
my @p_sol=();
my @u_sol=();


&clean_start;
foreach $e (@EIn) {
    foreach $rho (@rhoIn) {
	&clean_dir;
	&set_sim_params($e,$rho);
	&make_dek($dekfile,$parfile);
     
	$status = system "nohup", "mpirun", "-np", "2", "flash3", "-par_file", $parfile, "&";
	die if ($status);
  
	$status = system $idl, "idlrun";
	die if ($status);

	&process_data($max_dumps);
    }  
}



sub clean_start {
  my @dlist = glob "dumps/*.log nohup.out dumps/aerror.dat dumps/sb_* dumps/rho_* dumps/p_* dumps/u_* dumps/rho-* dumps/p-* dumps/u-*";
  foreach (@dlist) {
    unlink $_;
  }
}


sub clean_dir {
  my @dlist = glob "dumps/sb_* dumps/rho_* dumps/p_* dumps/u_* dumps/rho-* dumps/p-* dumps/u-*";
  foreach (@dlist) {
    unlink $_;
  }
}


sub set_sim_params {
  my $eIn    = shift;
  my $rhoIn  = shift;

  $sim_param{"sim_EIn"} = $eIn;
  $sim_param{"sim_rhoIn"} = $rhoIn;
}


sub make_dek {
  my $sfilename = shift;
  my $dfilename = shift;

  open INDEK, "<$sfilename" or die;
  open OUTDEK, ">$dfilename" or die;

  while (<INDEK>) {
      chomp;

      if (/sim_EIn/i) {
	  print OUTDEK "sim_EIn = ", $sim_param{"sim_EIn"} . "\n";
      }
      elsif (/sim_rhoIn/i) {
	  print OUTDEK "sim_rhoIn = ", $sim_param{"sim_rhoIn"} . "\n";
      }
      else {
	  print OUTDEK "$_\n";
      }
  }

  close INDEK;
  close OUTDEK;
}





sub collect_dumps {
    my $dnum = shift;
    my $n;
    my $t;
    my $i;

    @r=(); 
    @rho=();
    @p=(); 
    @u=();
    
    $i = 0;
    open RHOFILE, "<dumps/rho_$dnum" or die;
    while (<RHOFILE>) {
	chomp;
	if (!/\#/i) {
	    ($x,$a)  = split(" ",$_);
	    $r[$i]   = $x;
	    $rho[$i] = $a;
            $i++; 
	}
	else {
	    ($s,$n,$t) = split(" ",$_);
	}
    }
    close RHOFILE;

    $i = 0;
    open PFILE, "<dumps/p_$dnum" or die;
    while (<PFILE>) {
	chomp;
	if (!/\#/i) {
	    ($x,$a)  = split(" ",$_);
	    $p[$i]   = $a;
            $i++; 
	}
    }
    close PFILE;

    $i = 0;
    open UFILE, "<dumps/u_$dnum" or die;
    while (<UFILE>) {
	chomp;
	if (!/\#/i) {
	    ($x,$a)  = split(" ",$_);
	    $u[$i]   = $a;
            $i++; 
	}
    }
    close UFILE;

    open RFILE, ">dumps/r.dat" or die;
    print RFILE "# $n\n";
    for ($i=0; $i<$n; $i++) {
	print RFILE "$r[$i]\n";
    }
    close RFILE;

    ($n,$t);
}



sub make_sol_params {
    my $t        = shift;
    my $sfilename = shift;
    my $dfilename = shift;

    open SFILE, "<$sfilename" or die;
    open DFILE, ">$dfilename" or die;

    my $i=1;
    while (<SFILE>) {
	chomp;
	if ($i==5) {
	    print DFILE "$t \#t\n";
	}
	else {
	    print DFILE "$_\n";
	}
	$i++;
    }
   
    close SFILE;
    close DFILE;
}


sub read_sol {
    my $filename = shift;
 
    @rho_sol = ();
    @p_sol = ();
    @u_sol = ();
  
    my $i = 0;
    open SOLFILE, "<$filename" or die;
    while (<SOLFILE>) {
	chomp;
	if (!/\#/i) {
	    ($x,$v,$dens,$pres,$velx)  = split(" ",$_);
	    $rho_sol[$i] = $dens;
	    $p_sol[$i]   = $pres;
            $u_sol[$i]   = $velx;
            $i++; 
	}
    }
    close SOLFILE;
}



sub calc_errors {
    my $t = shift;
    my $filename = shift;

    my $rho_sum = 0.0;
    my $p_sum = 0.0;
    my $u_sum = 0.0;

    my $srho = 0.0;
    my $sp = 0.0;
    my $su = 0.0;
    
    my $i;
    my $n;

    my $r1;
    my $r2;
    my $dr;

    $n = @r;
    for ($i=0; $i<($n-1); $i++) {
	if ($i==0) {
	    $r1 = 0.0;
	}
	$r2 = 0.5*($r[$i+1] + $r[$i]);
	$dr = $r2 - $r1;

	$rho_sum = $rho_sum + abs($rho[$i] - $rho_sol[$i])*$dr;
	$p_sum = $p_sum + abs($p[$i] - $p_sol[$i])*$dr;
	$u_sum = $u_sum + abs($u[$i] - $u_sol[$i])*$dr;

	$srho = $srho + $rho_sol[$i]*$dr;
	$sp = $sp +  $p_sol[$i]*$dr;
	$su = $su + $u_sol[$i]*$dr;

	$r1 = $r2;
    }
    
    $rho_sum = $rho_sum / $srho;
    $p_sum = $p_sum / $sp;
    $u_sum = $u_sum / $su;

    open EFILE, ">>$filename" or die;
    print EFILE "$t $rho_sum $p_sum $u_sum\n";
    close EFILE;
}



sub plot_vars {
    my $dnum = shift;
    my $solfile = shift;

    my $zdnum = $dnum;
    if ($zdnum < 1000) {
	$zdnum = 0 . $zdnum;
    }
    if ($zdnum < 100) {
	$zdnum = 0 . $zdnum;
    }
    if ($zdnum < 10) {
	$zdnum = 0 . $zdnum;
    }

    open PFILE, ">>dumps/plot_data" or die;

    print PFILE "set term $term\n";
    print PFILE "set key top right\n";
    print PFILE "set key box\n";
    print PFILE "set output \"dumps/rho_$zdnum" . $dotterm . "\"\n";
    print PFILE "set title \"$title: Density\"\n";
    print PFILE "set xlabel \"r\"\n";
    print PFILE "set ylabel \"rho\"\n";
    print PFILE "plot \"$solfile\" using 1:3 w l lw $lw lt -1 t \"Solution\", ";
    print PFILE "\"dumps/rho_$dnum\" using 1:2 w l lw $lw lt 1 t \"Simulation\"\n";

    print PFILE "set term $term\n";
    print PFILE "set key top right\n";
    print PFILE "set key box\n";
    print PFILE "set output \"dumps/p_$zdnum" . $dotterm . "\"\n";
    print PFILE "set title \"$title: Pressure\"\n";
    print PFILE "set xlabel \"r\"\n";
    print PFILE "set ylabel \"p\"\n";
    print PFILE "plot \"$solfile\" using 1:4 w l lw $lw lt -1 t \"Solution\", ";
    print PFILE "\"dumps/p_$dnum\" using 1:2 w l lw $lw lt 1 t \"Simulation\"\n";

    print PFILE "set term $term\n";
    print PFILE "set key top right\n";
    print PFILE "set key box\n";
    print PFILE "set output \"dumps/u_$zdnum" . $dotterm . "\"\n";
    print PFILE "set title \"$title: Velocity\"\n";
    print PFILE "set xlabel \"r\"\n";
    print PFILE "set ylabel \"u\"\n";
    print PFILE "plot \"$solfile\" using 1:5 w l lw $lw lt -1 t \"Solution\", ";
    print PFILE "\"dumps/u_$dnum\" using 1:2 w l lw $lw lt 1 t \"Simulation\"\n";

    close PFILE;
}

sub process_data {
    my $max_dnum = shift;

    my $sedovparam = "dumps/sedov.param";
    my $paramfile  = "dumps/this.param";
    my $errfile    = "dumps/aerror.dat";

    my $i;

    for ($i=0; $i<=$max_dnum; $i++) {
	($n,$t) = &collect_dumps($i);
	
	&make_sol_params($t,$sedovparam,$paramfile);

	$status = system $sedovsolver, "$paramfile", "-v", "-alpha=0.851", "-rfile=dumps/r.dat";
	die if ($status);

	my $solfile  = "dumps/sedov_" . $i . ".dat";	
	system "mv", "sedov.dat", $solfile;

	&read_sol($solfile);
	&calc_errors($t,$errfile);
	&plot_vars($i,$solfile);
    }

    system "gnuplot", "dumps/plot_data";
}
