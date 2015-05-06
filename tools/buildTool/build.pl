#!/usr/bin/perl -w
BEGIN {
  (my $tmp = $0) =~ s/[^\/]+$//g;
  push(@INC,$tmp);
}
use config;
   @flist = ();
my @strippedFList = ();
my @objlist = ();
my $reflist;
my $depref;
my $mioref;
my $objref;

my %deplist = ();
my %modsInObj = ();
my %objOfMod = ();

printf "Starting Fortran build process...\n";
foreach my $sourceFile (@config::sourcefiles) { require $sourceFile; }
mkdir($config::buildPath);
foreach my $file (@flist) {
  (my $outFile = $file) =~ s/\//_/g;
  `$config::PRECOMP $config::PRECOMPFLAGS $config::path/$file -o $config::buildPath/$outFile`;
  `touch -r $config::path/$file $config::buildPath/$outFile`;
  push(@strippedFList,$outFile); 
}

printf "Copying Skript Files...\n";
`cp $config::preScript $config::buildPath/`;
`cp $config::postScript $config::buildPath/`;
(my $tmp = $config::postScript) =~ s/.*\///;
`perl -pi -e 's/PLACEHOLDER_COMPILER_ID/$config::CompilerID/g' $config::buildPath/$tmp`;
`cp $config::compModFileScript $config::buildPath/`;
($depref, $mioref, $objref) = buildModuleHash(\@strippedFList, $config::buildPath, \%deplist, \%modsInObj, \%objOfMod);
%deplist = %$depref;
%modsInObj = %$mioref;
%objOfMod = %$objref;

(my $preScript  = $config::preScript ) =~ s/.*\///;
(my $postScript = $config::postScript) =~ s/.*\///;

printf "Writing Makefile...\n";
open(MAKEFILE,">$config::buildPath/makefile");
printf MAKEFILE "FC=$config::FC\n";
printf MAKEFILE "FCFLAGS=$config::FCFLAGS\n";
printf MAKEFILE "LD=$config::LD\n";
printf MAKEFILE "INCLUDES=$config::INCLUDES\n";
printf MAKEFILE "LIBS=$config::LIBS\n";
printf MAKEFILE "all: $config::programName\n";
foreach (@strippedFList) {
 my $fname = $_;
 my $oname = $fname;
 $oname =~ s/\.[fF]90$/\.o/;
 push(@objlist, $oname);
 my @deps = @{ $deplist{$fname} };
 my @mods = @{ $modsInObj{$fname} };
 my @mfiles = ();
 my %seen = ();
 printf MAKEFILE "$oname";
 foreach (@mods) {
   my $mf = modfile($_, $objOfMod{$_});
   push(@mfiles, $mf) unless $seen{$mf}++;
   printf MAKEFILE " $mf" unless ($mf eq $oname);
 }
 printf MAKEFILE ": $fname";
 foreach (@deps) {
   my $depend = $_;
   my $dm = modfile($depend, $objOfMod{$depend});
   printf MAKEFILE " \\\n     $dm";
 }
 printf MAKEFILE "\n";
 printf MAKEFILE "\t\@sh $preScript @mfiles\n";
 printf MAKEFILE "\t\$(FC) $config::FCFLAGS $config::INCLUDES -c $fname\n";
 printf MAKEFILE "\t\@sh $postScript @mfiles\n";
 printf MAKEFILE "\n";
}

printf MAKEFILE "OBJS =";
foreach (@objlist)
{
    printf MAKEFILE " \\\n     $_";
}
printf MAKEFILE "\n";

printf MAKEFILE "$config::programName: \$(OBJS)\n";
printf MAKEFILE "\t\$(LD) $config::FCFLAGS $config::INCLUDES \$(OBJS) -o $config::programName $config::LIBS\n";
printf MAKEFILE "\n";

close(MAKEFILE);

printf "Compiling Program $config::programName...\n";
chdir($config::buildPath);
system("make");

printf "Finished!\n";

sub modfile
{
    # parameters
    # call with "$modname = modfile($mod, $objOfMod);"
    my ($mod, $oom) = @_;
    #
    #

    my $mfname = $mod.'.mod';
    return $mfname;
}

sub buildModuleHash
{
    # parameters
    # call with "($depref, $mioref, $objref) = buildModuleHash(\@flist, $config::buildPath, \%deplist, \%modsInObj, \%objOfMod);"
    my ($fref, $path, $depref, $mioref, $objref) = @_;
    #
    #

    my %deplist = %{$depref};

    my %objOfMod = %{$objref};

    my %modsInObj = %{$mioref};

    my @flist = @{$fref};

    foreach (@flist) {
	my $fname = $_;
	my $oname = $fname;
	$oname =~ s/\.f90$/\.o/;
	open(my $in,  "<",  "$path/$fname") 
	    or die "Can't open $path/$fname: $!";
	my @lines = <$in>;
	close $in or die "$in: $!";
	my @usedMods = ();
        my @modList = ();
	my $mc=0;
        my $nmod=0;
	foreach my $i (0 .. $#lines) {
	    if ($lines[$i] =~ /^ *module /i){
		my $modname = lc($lines[$i]);
		$modname =~ s/^ *module *(\S*).*$/$1/i;
		chomp $modname;
		if ( $modname ne "procedure" ) {
		    $objOfMod{$modname} = $oname;
                    $modList[$nmod] = $modname;
                    $nmod++;
		}
	    } elsif ($lines[$i] =~ /^ *use /i) {
		$usedMods[$mc] = lc($lines[$i]);
		$usedMods[$mc] =~ s/^ *use ([^\s,;]*).*$/$1/i;
		chomp($usedMods[$mc]);
		$mc++;
	    }
	}
	my %seen = ();
	my @uniq = ();
	foreach my $item (@usedMods) {
	    push(@uniq, $item) unless $seen{$item}++;
	}

	$deplist{$fname} = [ @uniq ];
        $modsInObj{$fname} = [ @modList ];

    }
    %{$depref} = %deplist;
    %{$objref} = %objOfMod;
    %{$mioref} = %modsInObj;

    return ($depref, $mioref, $objref);

}
