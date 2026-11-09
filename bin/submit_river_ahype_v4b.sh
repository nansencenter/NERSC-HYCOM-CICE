#!/bin/bash
#SBATCH --job-name=riverh_ahype_v4b
#SBATCH --account=nn2993k
#SBATCH --time=01:00:00
#SBATCH --qos=preproc
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH -o log/riverh_ahype_v4b.%j.out

# AHYPE v4b river forcing riverh.[ab] as a Slurm job. Gives the same file as
# river_ahype_v4b.sh. The spreading weights are computed once and reused for every year,
# so the whole period runs in one job (no job array).
#
# Usage, from the experiment directory. Run it directly, not with sbatch: it submits the
# Slurm job itself (the #SBATCH lines above are for that job):
#   $HOME/NERSC-HYCOM-CICE/bin/submit_river_ahype_v4b.sh [-s spreading] [first_year] [last_year]
# -s: spreading (bin/river_spreading.py), pointsource | coastal | diffuse (default pointsource).

spreading=pointsource
while getopts "s:h" opt ; do
    case $opt in
        s) spreading=$OPTARG ;;
        *) echo "Usage (from the experiment directory): $(basename $0) [-s spreading] [first_year] [last_year]"
           echo "   spreading: pointsource | coastal | diffuse (default pointsource); years default to 2002 2022"
           exit 1 ;;
    esac
done
shift $((OPTIND-1))
first=${1:-2002}
last=${2:-2022}
MODE=${3:-submit}   # submit (user) or run (the Slurm job submitted by this script)
if [ ! -f EXPT.src ]; then
    echo "Could not find EXPT.src. This script must be run in the experiment directory"
    exit 1
fi

if [ "$MODE" == "submit" ]; then
    # Inside a Slurm job $0 is a copy of this script, so pass the real bin directory on
    export AHYPE_BINDIR=$(cd $(dirname $0) && pwd)
    mkdir -p log
    jid=$(sbatch --parsable $AHYPE_BINDIR/$(basename $0) -s $spreading $first $last run) || exit 1
    echo "Job $jid: AHYPE v4b, spreading=$spreading, $first-$last -> force/rivers/<E>/riverh.[ab]"
    echo "Log: log/riverh_ahype_v4b.$jid.out"
    exit 0
fi

set --    # conda's activate would otherwise see the positional arguments
source ${HOME}/NERSC-HYCOM-CICE/environment/betzy_env.sh
module load Miniforge3/24.1.2-0
source ${EBROOTMINIFORGE3}/bin/activate
conda activate hycom-cice
$AHYPE_BINDIR/river_ahype_v4b.sh -s $spreading $first $last
