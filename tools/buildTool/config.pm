package config;

use strict;

$config::FC                = "ifort";
$config::FCFLAGS           = "-r8 -i4 -traceback -pg -vec-report0 -assume bscc -check bounds -check uninit -g";
$config::PRECOMP           = "$config::FC";
$config::PRECOMPFLAGS      = "-DMAXWELL3D -DPARTICLES -DRK -fpp -P";
$config::INCLUDES          = "-I ../CGNS";
$config::LIBS              = "-L ../CGNS/LINUX64 -lcgns";

$config::LD                = "$config::FC";

$config::CompilerID        = "INTEL-ifort9-on-LINUX";

@config::sourcefiles       = ("Particles.pm");

$config::path              = "./";
$config::buildPath         = "build";

$config::pwd               = `pwd`;
$config::command           = $0;
($config::skriptPath = $0) =~ s/[^\/]+$//g;
$config::preScript         = $config::skriptPath."prehandle_mods.sh";
$config::postScript        = $config::skriptPath."posthandle_mods.sh";
$config::compModFileScript = $config::skriptPath."compare_module_file.pl";

$config::programName       = "HALO";

1;
