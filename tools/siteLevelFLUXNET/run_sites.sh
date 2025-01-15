# First, we get a map of the repository by establishing where this script is located, then
# deducing where the root of the repository is.
script_location="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
rootdir=${script_location%%/tools*}
container=$rootdir/CLASSIC_container.simg
# Iterate through the sites and start a run for each one sequentially (using the container)
rm  timing.log
for f in $rootdir/inputFiles/FLUXNETsites/*; do
  if [ -d $f ]; then
    current=${f##*/}
    echo
    echo
    echo "Running $current..."
    echo
    #singularity exec $container $rootdir/bin/CLASSIC_serial $f/job_options_file.txt 0/0/0/0
    start=`date +%s.%N`
    $rootdir/bin/CLASSIC_serial $f/job_options_file.txt 0/0/0/0
    end=`date +%s.%N`
    runtime=$( echo "$end - $start" | bc -l )
    echo "Runtime for $current = $runtime"
    printf "Runtime for $current = $runtime\n" >> timing.log

  fi
done
