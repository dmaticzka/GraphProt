package StructureLibrary::RNAtoolsConfig;
#require Exporter;
#@ISA = qw(Exporter);
#@EXPORT = qw(
#			);
#@EXPORT_OK = qw(
#);

use strict;
use vars qw(%CONFIG);

# you may need to change some of these variables to fit to your system
our %CONFIG = (
    # All the configurations go here
    VERSION => '0.01',
    TMPDIR => 'tmp/',
    LIB_DIR => '/home/sita/Projects/Eclipse_Workspace/RNAtools/StructureLibrary/',
    RNAFOLD => '/usr/local/vrna/1.8.4/bin/RNAfold',
    RNAFOLDPARAMS => ' -p -d2 -noLP ',
    RNASHAPES => '/usr/local/rnashapes/2.1.1/bin/RNAshapes',
    RT	=>	0.61632, # RT = 1.98717 cal/mol/K * 310.15 = 0.61632 kcal/mol
    AUTHORS  => ["Sita Lange"],
    STR_FORMATS   => {
                   DOTBRACKET => '().-',
                   CONSTRAINT => '().|x<>',
                  },
);
1;