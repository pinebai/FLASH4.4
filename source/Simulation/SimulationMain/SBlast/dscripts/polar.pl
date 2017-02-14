#!/usr/bin/perl -w

my $title    = "Polar Plot";
# my $term     = "postscript eps color";
# my $dotterm  = ".eps";

my $term     = "png";
my $dotterm  = ".png";

my $visitbin = "/home/ahaque/visit/bin/visit"; # Path to visit
my $sedovsolver = "./sedov";                   # Path to sedov analytical solver

my @r=(); 
my @rho=();
my @p=(); 
my @u=();
my @ut=();
my @ux=();
my @uy=();

my @rho_sol=();
my @p_sol=();
my @u_sol=();

my $xmax = 1;
my $ymax = 1;

my $rmin_check = 0.01;
my $pmin_check = 0.01;

my $lw       = 2;
my @lt       = (5,4,3,2,1);

my $n_times  = 100;
my $n_angles = 36;

my $prefix = "sbpol";
my $sprefix = "sb";

my @mindex = (20,40,60,80,100);
my @subset = (0,1,2);
my %tlist;

my @adata;
my @rdata;

my $pi = 3.14159265358979;
my $dtheta =  $pi / $n_angles;



my $doall=0;

@dlist = glob "error* *.png";
foreach (@dlist) {
    unlink $_;
}


die if system $visitbin, "-default_format", "FLASH", "-cli", "-nowin", "-s", "visit_polar_dens.py";
die if system $visitbin, "-default_format", "FLASH", "-cli", "-nowin", "-s", "visit_polar_pres.py";
die if system $visitbin, "-default_format", "FLASH", "-cli", "-nowin", "-s", "visit_polar_velx.py";
die if system $visitbin, "-default_format", "FLASH", "-cli", "-nowin", "-s", "visit_polar_vely.py";

if (!$doall) {
foreach $it (@mindex) {
    my $filename = $prefix . "_pres_" . $it . "_0.out";
    $outname = $prefix . "_" . $it . ".dat";
    open OUTFILE, ">$outname" or die;

    ($t,$theta,$r) = &find_r($filename);
    $tlist{$it} = $t; 
    print OUTFILE "# $t\n";
    print OUTFILE "$theta $r\n";
    $adata[0] = $theta;
    $rdata[0] = $r;

    for ($ia=1; $ia<=$n_angles; $ia++) {
	$filename = $prefix . "_pres_" . $it . "_" . $ia . ".out";
	($t,$theta,$r) = &find_r($filename);
	print OUTFILE "$theta $r\n";
	$adata[$ia] = 3*$pi/2 - $ia*$dtheta;
	$rdata[$ia] = $r;
    }
    for ($ia=($n_angles-1); $ia>=0; $ia--) {
	print OUTFILE "$adata[$ia] $rdata[$ia]\n";
    }
    close OUTFILE;
    $plotname = $prefix . "_" . $it . "$dotterm";
    &make_polar($outname,$plotname);

    $outname = $prefix . "_r_" . $it . ".dat";
    open OUTFILE, ">$outname" or die;
    my $efilename = "error_" . $it . ".dat";    
    for ($ia=0; $ia<=$n_angles; $ia++) {
	my $pres_filename = $prefix . "_pres_" . $it . "_" . $ia . ".out";
        $t = read_pres_data($pres_filename);
	($t,$theta,$r) = &find_r($pres_filename);
	my $thetadeg = &rad2deg($theta);
	$r_err = $r - &r_sol($t);
	print OUTFILE "$thetadeg $r $r_err\n";

	my $dens_filename = $prefix . "_dens_" . $it . "_" . $ia . ".out";
        $t = read_dens_data($dens_filename);

	my $velx_filename = $prefix . "_velx_" . $it . "_" . $ia . ".out";
        $t = read_velx_data($velx_filename);

	my $vely_filename = $prefix . "_vely_" . $it . "_" . $ia . ".out";
        $t = read_vely_data($vely_filename);	

	&make_u_data($theta,"velr.dat");

	&make_sol_params($t,"sedov.param","this.param");
	&make_r_file("r.dat");

	$status = system $sedovsolver, "this.param", "-v", "-alpha=0.851", "-rfile=r.dat";
	die if ($status);
	&read_sol("sedov.dat");

	&make_plots($t,$thetadeg,"sedov.dat",$dens_filename,$pres_filename,"velr.dat",$it . "_" . $ia);
	&calc_errors($t,$theta,$efilename);
    }
    close OUTFILE;
}

&make_replots("rerror");
&make_eplots("error");
&make_mpolar($prefix . "_mpolar" . $doterm);
}
die("done for now");



#####################################################################################################
# Still in development

if ($doall) {
die if system $visitbin, "-default_format", "FLASH", "-cli", "-nowin", "-s", "visit_subset_pres.py";
die if system $visitbin, "-default_format", "FLASH", "-cli", "-nowin", "-s", "visit_subset_dens.py";
die if system $visitbin, "-default_format", "FLASH", "-cli", "-nowin", "-s", "visit_subset_velx.py";
die if system $visitbin, "-default_format", "FLASH", "-cli", "-nowin", "-s", "visit_subset_vely.py";

foreach $ia (@subset) {
    $filename = $sprefix . "_pres_0_" . $ia . ".out";
    $outname = $sprefix . "_rt_" . $ia . ".dat";
    open OUTFILE, ">$outname" or die;

    ($t,$theta,$r) = &find_r($filename);
    print OUTFILE "# Angle $theta\n";
    print OUTFILE "$t $r\n";   

    for ($it=0; $it<=$n_times; $it++) {
	$filename = $sprefix . "_pres_" . $it . "_" . $ia . ".out";
	($t,$theta,$r) = &find_r($filename);
	print OUTFILE "$t $r\n";

	my $pres_filename = $prefix . "_pres_" . $it . "_" . $ia . ".out";
        $t = read_pres_data($pres_filename);
	($t,$theta,$r) = &find_r($pres_filename);
	my $thetadeg = &rad2deg($theta);
	$r_err = $r - &r_sol($t);
	print OUTFILE "$thetadeg $r $r_err\n";

	my $dens_filename = $prefix . "_dens_" . $it . "_" . $ia . ".out";
        $t = read_dens_data($dens_filename);

	my $velx_filename = $prefix . "_velx_" . $it . "_" . $ia . ".out";
        $t = read_velx_data($velx_filename);

	my $vely_filename = $prefix . "_vely_" . $it . "_" . $ia . ".out";
        $t = read_vely_data($vely_filename);	

	&make_u_data($theta,"velr.dat");

	&make_sol_params($t,"sedov.param","this.param");
	&make_r_file("r.dat");

	$status = system $sedovsolver, "this.param", "-v", "-alpha=0.851", "-rfile=r.dat";
	die if ($status);
	&read_sol("sedov.dat");

	&make_plots($t,$thetadeg,"sedov.dat",$dens_filename,$pres_filename,"velr.dat",$it . "_" . $ia);
    }

    close OUTFILE;
}
}

###################################################################################


sub rad2deg {
    my $rad = shift;
    my $angle;

    $angle = 180.0 * $rad / $pi;

    $angle = sprintf("%.0f",$angle);
}


sub r_sol {
    my $t = shift;

    my $E0 = 1;
    my $alpha = 0.851;
    my $rho = 1.0;
    
    my $r = ($E0/($rho*$alpha))**0.2 * $t**0.4;

    $r;
}

sub read_dens_data {
    my $filename = shift;
    my $t;
    my $i;
    my $line;

    @r=();
    @rho=();

    $line = 0;
    $i = 0;
    open RHOFILE, "<$filename" or die;
    while (<RHOFILE>) {
	chomp;
	if (!/\#/i) {
	    ($x,$a)  = split(" ",$_);
	    $r[$i]   = $x;
	    $rho[$i] = $a;
            $i++; 
	}
	else {
	    if ($line==0) {
		($s,$t) = split(" ",$_);
	    }
	}
	$line++;
    }
    close RHOFILE;

    return $t;
}


sub read_pres_data {
    my $filename = shift;
    my $t;
    my $i;
    my $line;

    @r=();
    @p=();

    $line = 0;
    $i = 0;
    open PFILE, "<$filename" or die;
    while (<PFILE>) {
	chomp;
	if (!/\#/i) {
	    ($x,$a)  = split(" ",$_);
	    $r[$i]   = $x;
	    $p[$i]   = $a;
            $i++; 
	}
	else {
	    if ($line==0) {
		($s,$t) = split(" ",$_);
	    }
	}
	$line++;
    }
    close PFILE;

    return $t;
}


sub read_velx_data {
    my $filename = shift;
    my $t;
    my $i;
    my $line;

    @r=();
    @ux=();

    $line = 0;
    $i = 0;
    open UXFILE, "<$filename" or die;
    while (<UXFILE>) {
	chomp;
	if (!/\#/i) {
	    ($x,$a)  = split(" ",$_);
	    $r[$i]   = $x;
	    $ux[$i]  = $a;
            $i++; 
	}
	else {
	    if ($line==0) {
		($s,$t) = split(" ",$_);
	    }
	}
	$line++;
    }
    close UXFILE;

    return $t;
}



sub read_vely_data {
    my $filename = shift;
    my $t;
    my $i;
    my $line;

    @r=();
    @uy=();

    $line = 0;
    $i = 0;
    open UYFILE, "<$filename" or die;
    while (<UYFILE>) {
	chomp;
	if (!/\#/i) {
	    ($x,$a)  = split(" ",$_);
	    $r[$i]   = $x;
	    $uy[$i]  = $a;
            $i++; 
	}
	else {
	    if ($line==0) {
		($s,$t) = split(" ",$_);
	    }
	}
	$line++;
    }
    close UYFILE;

    return $t;
}



sub make_u_data {
    my $theta = shift;
    my $filename = shift;
    my $n;
    my $i;

    $n = @ux;

    my $n_x = cos($theta);
    my $n_y = sin($theta);
    my $t_x = cos($theta + $pi/2);
    my $t_y = sin($theta + $pi/2);

    open UFILE, ">$filename" or die;
    for ($i=0; $i<$n; $i++) {

	$u[$i]  = $ux[$i]*$n_x + $uy[$i]*$n_y;
	$ut[$i] = $ux[$i]*$t_x + $uy[$i]*$t_y;

	print UFILE "$r[$i] $u[$i] $ut[$i]\n";
    }
    close UFILE;
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


sub make_r_file {
    my $filename = shift;
    my $n;
    my $i;

    $n = @r;
    open RFILE, ">$filename" or die;
    print RFILE "# $n\n";
    for ($i=0; $i<$n; $i++) {
	print RFILE $r[$i] . "\n";
    }
    close RFILE;
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
    my $angle = shift;
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

    $angle = &rad2deg($angle);
    open EFILE, ">>$filename" or die;
    print EFILE "$t $angle $rho_sum $p_sum $u_sum\n";
    close EFILE;
}



sub make_plots {
    my $t = shift;
    my $a = shift;
    my $solfile = shift;
    my $dfile = shift;
    my $pfile = shift;
    my $ufile = shift;
    my $suffix = shift;

    my $ofile;

    $t = sprintf("%4.2f",$t);
    $a = sprintf("%5.2f",$a);

    open PLOTFILE, ">ta_plotfile" or die;
    
    $ofile = "dens_" . $suffix . $dotterm;
    print PLOTFILE "set term $term\n";
    print PLOTFILE "set output \"$ofile\"\n";
    print PLOTFILE "set title \"Sedov: Density, t=$t, angle=$a\"\n";
    print PLOTFILE "set xlabel \"r\"\n";
    print PLOTFILE "set ylabel \"rho\"\n";
    print PLOTFILE "set key top right\n";
    print PLOTFILE "set key box\n";
    print PLOTFILE "plot \"$solfile\" using 1:3 w l lt -1 lw $lw t \"Solution\", ";
    print PLOTFILE "\"$dfile\" using 1:2 w l lt 1 lw $lw t \"Simulation\"\n";

    $ofile = "pres_" . $suffix . $dotterm;
    print PLOTFILE "set term $term\n";
    print PLOTFILE "set output \"$ofile\"\n";
    print PLOTFILE "set title \"Sedov: Pressure, t=$t, angle=$a\"\n";
    print PLOTFILE "set xlabel \"r\"\n";
    print PLOTFILE "set ylabel \"p\"\n";
    print PLOTFILE "set key top right\n";
    print PLOTFILE "set key box\n";
    print PLOTFILE "plot \"$solfile\" using 1:4 w l lt -1 lw $lw t \"Solution\", ";
    print PLOTFILE "\"$pfile\" using 1:2 w l lt 1 lw $lw t \"Simulation\"\n";

    $ofile = "velr_" . $suffix . $dotterm;
    print PLOTFILE "set term $term\n";
    print PLOTFILE "set output \"$ofile\"\n";
    print PLOTFILE "set title \"Sedov: Radial Velocity, t=$t, angle=$a\"\n";
    print PLOTFILE "set xlabel \"r\"\n";
    print PLOTFILE "set ylabel \"u_r\"\n";
    print PLOTFILE "set key top right\n";
    print PLOTFILE "set key box\n";
    print PLOTFILE "plot \"$solfile\" using 1:5 w l lt -1 lw $lw t \"Solution\", ";
    print PLOTFILE "\"$ufile\" using 1:2 w l lt 1 lw $lw t \"Simulation\"\n";

    $ofile = "velt_" . $suffix . $dotterm;
    print PLOTFILE "set term $term\n";
    print PLOTFILE "set output \"$ofile\"\n";
    print PLOTFILE "set title \"Sedov: Tangential Velocity, t=$t, angle=$a\"\n";
    print PLOTFILE "set xlabel \"r\"\n";
    print PLOTFILE "set ylabel \"u_t\"\n";
    print PLOTFILE "set key top right\n";
    print PLOTFILE "set key box\n";
    print PLOTFILE "plot 0 w l lt -1 lw $lw t \"Solution\", ";
    print PLOTFILE "\"$ufile\" using 1:3 w l lt 1 lw $lw t \"Simulation\"\n";

    close PLOTFILE;

    my $status = system "gnuplot", "<", "ta_plotfile";
    die if $status;
}


sub make_replots {
    my $pfile = shift;
    my $plotfile = "re_plotfile";
    my $dfile;
    my $ofile;

    my $i;
    my $n = @mindex;

    $ofile = $pfile . $dotterm;
    open PLOTFILE, ">$plotfile" or die;
    print PLOTFILE "set term $term\n";
    print PLOTFILE "set output \"$ofile\"\n";
    print PLOTFILE "set title \"Sedov: Shock Location Error vs Angle\"\n";
    print PLOTFILE "set ylabel \"Shock Location Error\"\n";
    print PLOTFILE "set xlabel \"Angle\"\n";
    print PLOTFILE "set xrange [-90:90]\n";
    print PLOTFILE "set key bottom right\n";
    print PLOTFILE "set key box\n";

    for ($i=0; $i<$n; $i++) {
        $dfile = $prefix . "_r_" . $mindex[$i] . ".dat";
	$tit = "t=" . sprintf("%.2f",$tlist{$mindex[$i]});
	if ($i==0) {
	    print PLOTFILE "plot \"$dfile\" using 1:3 w lp lt $lt[$i] lw $lw t \"$tit\"";
	}
	else {
	    print PLOTFILE ", \"$dfile\" using 1:3 w lp lt $lt[$i] lw $lw t \"$tit\"";
	}
    }
    print PLOTFILE "\n";

    close PLOTFILE;

    my $status = system "gnuplot", "<", $plotfile;
    die if $status;
}


sub make_eplots {
    my $pfile = shift;
    my $plotfile = "plotfile";
    my $dfile;
    my $ofile;

    my $i;
    my $n = @mindex;

    $ofile = $pfile . "_dens" . $dotterm;
    open PLOTFILE, ">$plotfile" or die;
    print PLOTFILE "set term $term\n";
    print PLOTFILE "set output \"$ofile\"\n";
    print PLOTFILE "set title \"Sedov: Density Error vs Angle\"\n";
    print PLOTFILE "set ylabel \"Scaled L1 Error\"\n";
    print PLOTFILE "set xlabel \"Angle\"\n";
    print PLOTFILE "set xrange [-90:90]\n";
    print PLOTFILE "set key bottom right\n";
    print PLOTFILE "set key box\n";

    for ($i=0; $i<$n; $i++) {
	$dfile = "error_" . $mindex[$i] . ".dat";
	$tit = "t=" . sprintf("%.2f",$tlist{$mindex[$i]});
	if ($i==0) {
	    print PLOTFILE "plot \"$dfile\" using 2:3 w lp lt $lt[$i] lw $lw t \"$tit\"";
	}
	else {
	    print PLOTFILE ", \"$dfile\" using 2:3 w lp lt $lt[$i] lw $lw t \"$tit\"";
	}
    }
    print PLOTFILE "\n";

    $ofile = $pfile . "_pres" . $dotterm;
    print PLOTFILE "set term $term\n";
    print PLOTFILE "set output \"$ofile\"\n";
    print PLOTFILE "set title \"Sedov: Pressure Error vs Angle\"\n";
    print PLOTFILE "set ylabel \"Scaled L1 Error\"\n";
    print PLOTFILE "set xlabel \"Angle\"\n";
    print PLOTFILE "set xrange [-90:90]\n";
    print PLOTFILE "set key bottom right\n";
    print PLOTFILE "set key box\n";
   
    for ($i=0; $i<$n; $i++) {
	$dfile = "error_" . $mindex[$i] . ".dat";
	$tit = "t=" .  sprintf("%.2f",$tlist{$mindex[$i]});
	if ($i==0) {
	    print PLOTFILE "plot \"$dfile\" using 2:4 w lp lt $lt[$i] lw $lw t \"$tit\"";
	}
	else {
	    print PLOTFILE ", \"$dfile\" using 2:4 w lp lt $lt[$i] lw $lw t \"$tit\"";
	}
    }
    print PLOTFILE "\n";

    $ofile = $pfile . "_velr" . $dotterm;
    print PLOTFILE "set term $term\n";
    print PLOTFILE "set output \"$ofile\"\n";
    print PLOTFILE "set title \"Sedov: Velocity Error vs Angle\"\n";
    print PLOTFILE "set ylabel \"Scaled L1 Error\"\n";
    print PLOTFILE "set xlabel \"Angle\"\n";
    print PLOTFILE "set xrange [-90:90]\n";
    print PLOTFILE "set key bottom right\n";
    print PLOTFILE "set key box\n";   

    for ($i=0; $i<$n; $i++) {
	$dfile = "error_" . $mindex[$i] . ".dat";
	$tit = "t=" .  sprintf("%.2f",$tlist{$mindex[$i]});
	if ($i==0) {
	    print PLOTFILE "plot \"$dfile\" using 2:5 w lp lt $lt[$i] lw $lw t\"$tit\"";
	}
	else {
	    print PLOTFILE ", \"$dfile\" using 2:5 w lp lt $lt[$i] lw $lw t \"$tit\"";
	}
    }
    print PLOTFILE "\n";

    close PLOTFILE;

    $status = system "gnuplot", "<", "$plotfile";
    die if $status;
}


sub find_r {
    my $filename = shift;

    my $t;
    my $theta;
    my $r_max = 0.0;
    
    print "--->$filename\n";

    open INFILE, "<$filename" or die;

    my $last_r = $r_max;
    my $looking = 1;

    while (<INFILE>) {
	chomp;
	
	if (/Angle/i) {
	    ($s1,$s2,$theta) = split(' ',$_);
	}
	elsif (/Flash/i) {
	    
	}
	elsif (/\#/i) {
	    ($s1,$t) = split(' ',$_);
	}
	else {
	    ($r,$p) = split(' ',$_);
	    if (($p < $pmin_check) && ($r > $rmin_check) && $looking) {
		$r_max = $last_r;
		$looking = 0;
	    }
	    $last_r = $r;
	}
    }

    close INFILE;

    ($t,$theta,$r_max);
}



sub make_polar {
    my $dfile = shift;
    my $pfile = shift;

    my $plotfile = "plotfile";
    
    open PLOTFILE, ">$plotfile" or die;
    print PLOTFILE "set term $term\n";
    print PLOTFILE "set output \"$pfile\"\n";
    print PLOTFILE "set title \"$title\"\n";
    print PLOTFILE "set xlabel \"x\"\n";
    print PLOTFILE "set ylabel \"y\"\n";
    print PLOTFILE "set nokey\n";
    print PLOTFILE "set polar\n";
    print PLOTFILE "set grid polar\n";
    print PLOTFILE "set size square\n";
    print PLOTFILE "set xrange [" . -$xmax . ":" . $xmax . "]\n";
    print PLOTFILE "set yrange [" . -$ymax . ":" . $ymax . "]\n";
    print PLOTFILE "plot \"$dfile\" with lines lt 1 lw $lw\n";
    close PLOTFILE;

    system "gnuplot", "<", "$plotfile";
}



sub get_time {
    my $filename = shift;
    my $t;

    open INFILE, "<$filename" or die;
    
    while (<INFILE>) {
	chomp;

	if (/\#/) {
	    ($s,$t) = split(' ',$_);
	}
    }
    close INFILE;

    $t;
}


sub make_mpolar {
    my $pfile = shift;
    my $t;

    my $plotfile = "plotfile";
    
    open PLOTFILE, ">$plotfile" or die;
    print PLOTFILE "set term $term\n";
    print PLOTFILE "set output \"$pfile\"\n";
    print PLOTFILE "set title \"$title\"\n";
    print PLOTFILE "set xlabel \"x\"\n";
    print PLOTFILE "set ylabel \"y\"\n";
    print PLOTFILE "set key bottom right\n";
    print PLOTFILE "set key box\n";
    print PLOTFILE "set polar\n";
    print PLOTFILE "set grid polar\n";
    print PLOTFILE "set size square\n";
    print PLOTFILE "set xrange [" . -$xmax . ":" . $xmax . "]\n";
    print PLOTFILE "set yrange [" . -$ymax . ":" . $ymax . "]\n";

    $ilt = 0;
    foreach $index (@mindex) {
	$dfile = $prefix . "_" . $index . ".dat";
	$t = &get_time($dfile);
	$t = sprintf("%4.2f",$t);
        if ($ilt==0) {
	    print PLOTFILE "plot \"$dfile\" with lines lt $lt[$ilt] lw $lw t \"t=$t\"";
	}
	else {
	    print PLOTFILE ", \"$dfile\" with lines lt $lt[$ilt] lw $lw t \"t=$t\"";
	}
	$ilt++;
    }

    print PLOTFILE "\n";
    close PLOTFILE;

    system "gnuplot", "<", "$plotfile";
}
