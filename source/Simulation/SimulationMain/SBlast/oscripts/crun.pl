#!/usr/bin/perl -w

# Convergence study for the Sedov problem

my $title    = "Sedov";
my $term     = "png";
my $dotterm  = ".png";
my $lw       = 2;

my $dekfile  = "flash1d_uni.par";
my $parfile  = "this.par";

my $idl      = "/usr/local/rsi/idl_6.1/bin/idl";   # Path to idl
my $sedovsolver = "dumps/sedov";                   # Path to Sedov analytical solver

my $max_dumps = 100;                               # Last checkpoint file number

my @colors = (5,4,2,3,1);                          # gnuplot line styles

my @n_blocks = (50,100,150,200,250);               # Number of blocks

my %sim_params;

my @r=(); 
my @rho=();
my @p=(); 
my @u=();

my @rho_sol=();
my @p_sol=();
my @u_sol=();

my $nb;


&clean_start;
foreach $nb (@n_blocks) {
    &clean_dir;
    &set_sim_params($nb);
    &make_dek($dekfile,$parfile);
     
    $status = system "nohup", "mpirun", "-np", "2", "flash3", "-par_file", $parfile, "&";
    die if ($status);
  
    $status = system $idl, "idlrunlast";
    die if ($status);

    &process_data($nb,$max_dumps);
}
&make_plots;


sub clean_start {
    &clean_dir;

    my @dlist = glob "dumps/*.log nohup.out dumps/error.dat dumps/sb_* dumps/*.dat dumps/rho_* dumps/p_* dumps/u_* dumps/rho-* dumps/p-* dumps/u-*";
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
  my $nb = shift;

  $sim_param{"nblocks"} = $nb;
}


sub make_dek {
  my $sfilename = shift;
  my $dfilename = shift;

  open INDEK, "<$sfilename" or die;
  open OUTDEK, ">$dfilename" or die;

  while (<INDEK>) {
      chomp;

      if (/nblockx/i) {
	  print OUTDEK "nblockx = ", $sim_param{"nblocks"} . "\n";
      }
#      elsif (/nblocky/i) {
#	  print OUTDEK "nblocky = ", $sim_param{"nblocks"} . "\n";	  
#      }
#      elsif (/nblockz/i) {
#	  print OUTDEK "nblockz = ", $sim_param{"nblocks"} . "\n";
#      }
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

    print "dnum = $dnum\n";
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
    my $nb = shift;
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

#    print "u_sum = $u_sum \n";
#    print "su = $su \n";
    
    $rho_sum = $rho_sum / $srho;
    $p_sum = $p_sum / $sp;
    $u_sum = $u_sum / $su;

#    print "norm = $u_sum \n";

    open EFILE, ">>$filename" or die;
    print EFILE "$nb $rho_sum $p_sum $u_sum\n";
    close EFILE;
}



sub process_data {
    my $nb = shift;
    my $max_dnum = shift;

    my $sedovparam = "dumps/sedov.param";
    my $paramfile  = "dumps/this.param";
    my $errfile    = "dumps/error.dat";

    ($n,$t) = &collect_dumps($max_dnum);
	
    &make_sol_params($t,$sedovparam,$paramfile);

    $status = system $sedovsolver, "$paramfile", "-v", "-alpha=0.851", "-rfile=dumps/r.dat";
    die if ($status);

    my $solfile  = "dumps/sedov.dat";
    system "mv", "sedov.dat", $solfile;

    &read_sol($solfile);
    &calc_errors($nb,$errfile);

    system "cp", "dumps/rho_$max_dnum", "dumps/rhonb" . $nb;
    system "cp", "dumps/p_$max_dnum", "dumps/pnb" . $nb;
    system "cp", "dumps/u_$max_dnum", "dumps/unb" . $nb;    
}



sub make_plots {
    my $solfile  = "dumps/sedov.dat";

    open PLOTFILE, ">dumps/plot_data" or die;
    print PLOTFILE "set term $term\n";
    print PLOTFILE "set key top left\n";
    print PLOTFILE "set key box\n";
    print PLOTFILE "set output \"dumps/$title" . "_rho" . "$dotterm\"\n";
    print PLOTFILE "set title \"$title: Density\"\n";
    print PLOTFILE "set xlabel \"r\"\n";
    print PLOTFILE "set ylabel \"rho\"\n";
    print PLOTFILE "plot \"$solfile\" using 1:3 w l lw $lw lt -1 t \"Solution\"";
    for ($i=0; $i<@n_blocks; $i++) {
	$dfile = "dumps/rhonb" . $n_blocks[$i];
	print PLOTFILE ", \"$dfile\" using 1:2 w l lw $lw lt $colors[$i]  t \"$n_blocks[$i] blocks\"";
    }
    close PLOTFILE;    
    system "gnuplot", "dumps/plot_data";
    unlink "dumps/plot_data";

    open PLOTFILE, ">dumps/plot_data" or die;
    print PLOTFILE "set term $term\n";
    print PLOTFILE "set key top left\n";
    print PLOTFILE "set key box\n";
    print PLOTFILE "set output \"dumps/$title" . "_p" . "$dotterm\"\n";
    print PLOTFILE "set title \"$title: Pressure\"\n";
    print PLOTFILE "set xlabel \"r\"\n";
    print PLOTFILE "set ylabel \"p\"\n";
    print PLOTFILE "plot \"$solfile\" using 1:4 w l lw $lw lt -1 t \"Solution\"";
    for ($i=0; $i<@n_blocks; $i++) {
	$dfile = "dumps/pnb" . $n_blocks[$i];
	print PLOTFILE ", \"$dfile\" using 1:2 w l lw $lw lt $colors[$i]  t \"$n_blocks[$i] blocks\"";
    }
    close PLOTFILE;    
    system "gnuplot", "dumps/plot_data";
    unlink "dumps/plot_data";

   open PLOTFILE, ">dumps/plot_data" or die;
    print PLOTFILE "set term $term\n";
    print PLOTFILE "set key top left\n";
    print PLOTFILE "set key box\n";
    print PLOTFILE "set output \"dumps/$title" . "_u" . "$dotterm\"\n";
    print PLOTFILE "set title \"$title: Velocity\"\n";
    print PLOTFILE "set xlabel \"r\"\n";
    print PLOTFILE "set ylabel \"u\"\n";
    print PLOTFILE "plot \"$solfile\" using 1:5 w l lw $lw lt -1 t \"Solution\"";
    for ($i=0; $i<@n_blocks; $i++) {
	$dfile = "dumps/unb" . $n_blocks[$i];
	print PLOTFILE ", \"$dfile\" using 1:2 w l lw $lw lt $colors[$i]  t \"$n_blocks[$i] blocks\"";
    }
    close PLOTFILE;    
    system "gnuplot", "dumps/plot_data";
    unlink "dumps/plot_data";
}
