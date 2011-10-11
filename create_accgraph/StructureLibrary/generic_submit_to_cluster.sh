#!/bin/bash
#$ -cwd
#$ -e /home/sita/Projects/sge.err/err.$JOB_ID
#$ -o /home/sita/Projects/sge.err/out.$JOB_ID
#$ -l h_vmem=3G
#$ -l nice_cpu=1
#$ -m ea
#$ -soft -l scratch1_local=1

# -l h_cpu=2:0:0 
# -l h_rt=3:0:0 


######################################################################
## job script "generic_submit_to_cluster.sh"
## so this is a generic shell script that can call a certain perl script
## with the given parameters. 
## This is useful if you need to commit several similar jobs to the cluster
######################################################################



#### all paths that need to be predefined
#export PERL5LIB=$HOME/install/bin:$HOME/install/lib/perl:$PERL5LIB
export PERL5LIB=/usr/bin/perl5.10.0:$HOME/Install/lib/perl5:$HOME/Projects/Eclipse_Workspace/CVSROOT_HOME/Projects/Generalscripts:$HOME/Projects/Eclipse_Workspace/RNAtools:$PERL5LIB
#standard path variable for any programs that are used
export PATH=/usr/local/vrna/1.8.2/bin/:$HOME/bin/:$PATH

# all parameters for the shell script
working_dir=$1
perlscript=$2
params=$3

if [ -z "$SGE_TASK_ID" ]; then
	echo "No SGE environment found! Set SGE_TASK_ID to 0!"
	let SGE_TASK_ID=0
fi

cd $working_dir
# Here we use perl as set by the path variable to run the script. If
# called directly using ./$perlscript, bash looks at the first line
# of the script to determine the executable to use. This is usually set
# to #!/usr/bin/perl. On the sge cluster this is currently
# version v5.8.8, released in early 2006, which may not be what you want.
perl ./$perlscript $params -jobid $SGE_TASK_ID

