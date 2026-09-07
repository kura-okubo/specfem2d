#!/bin/bash
###################################################
#
# Validates the external-source coupling end to end:
#
#   1. mesher with COUPLING_IN = .true.  -> the coupling elements
#   2. a point-source run                -> the "true" wavefield
#   3. that wavefield, sampled at the coupling elements, becomes the injected
#      source
#   4. a run with the built-in source switched off, driven only by the injection
#   5. the two runs must agree at receivers outside the coupling surface
#
# There is no stored reference solution on purpose: the test compares the
# coupling against the point-source run computed in the same job, which is
# exactly the statement the method has to satisfy (see
# EXAMPLES/Validation/note/Note_of_Specfem2DCoupling.ipynb).
#
# EXAMPLES/Validation/FullSpace is the same setup, at full resolution.
#
###################################################

# example the DATA/ is taken from
NAME="Validation/FullSpace"

# relative location of SPECFEM2D EXAMPLES/ directory
EXAMPLES="../../EXAMPLES/"

# receivers used for the comparison: radius, and azimuths from north
R_VALIDATE=10000.0
AZIMUTHS="0 90 180 270"

# a broken coupling gives a correlation near zero. The threshold is loose on
# purpose: the injection leaves a high-frequency numerical oscillation that the
# note's own figures low-pass away before comparing, and that costs ~0.1 here.
MIN_CORRELATION=0.8

###################################################

testdir=`pwd`

# title
echo >> $testdir/results.log
echo "coupling validation ($NAME) in: $testdir" >> $testdir/results.log
echo >> $testdir/results.log

# checks if compilation done
if [ ! -e ./bin/xspecfem2D ]; then
  echo "no binaries found, compilation must be done before, please check..." >> $testdir/results.log
  exit 1
fi

# checks if example directory exists
if [ ! -e $EXAMPLES/$NAME/DATA/Par_file ]; then
  echo "example DATA in $EXAMPLES/$NAME not found, please check..." >> $testdir/results.log
  exit 1
fi

# how many MPI ranks: the injection is averaged across MPI interfaces, so run
# more than one rank when we can, to exercise that path too
NPROC=2
if [ ! -e ./bin/xmeshfem2D ]; then NPROC=1; fi
command -v mpirun > /dev/null 2>&1 || NPROC=1
# the build only speaks MPI when configure was given --with-mpi
grep -q "^COND_MPI_FCFLAGS *=.*WITH_MPI" ./Makefile 2>/dev/null || NPROC=1

# sed -i differs between GNU and BSD; upstream's tests assume GNU, this does not
sedi(){
  if sed --version > /dev/null 2>&1; then sed -i "$@"; else sed -i '' "$@"; fi
}

run_solver(){
  # $1 = executable
  if [ "$NPROC" -eq 1 ]; then
    ./bin/$1 >> $testdir/results.log 2>&1
  else
    mpirun -np $NPROC ./bin/$1 >> $testdir/results.log 2>&1
  fi
}

# cleanup
rm -rf ./DATA ./OUTPUT_FILES ./OUTPUT_FILES_P ./OUTPUT_FILES_C ./extsource
mkdir -p OUTPUT_FILES

# setup
cp -rp $EXAMPLES/$NAME/DATA .

###################################################
# shrink it so the whole chain runs in seconds:
# 50 km box at 500 m elements (Vs = 1732 m/s, source f0 = 0.2 Hz, so still far
# better than resolved), 14 s of propagation, receivers at 10 km
###################################################
sedi "s:^NPROC .*:NPROC                           = $NPROC:" DATA/Par_file
sedi "s:^NSTEP .*:NSTEP                           = 1400:" DATA/Par_file
sedi "s:^DT .*:DT                              = 1.0d-2:" DATA/Par_file
sedi "s:^nx .*:nx                              = 100:" DATA/Par_file
sedi "s:^1 200 1 200 1:1 100 1 100 1:" DATA/Par_file
sedi "s:^dR_ext .*:dR_ext                          = 500.0:" DATA/Par_file
sedi "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO      = 500:" DATA/Par_file
sedi "s:^NTSTEP_BETWEEN_OUTPUT_SEISMOS .*:NTSTEP_BETWEEN_OUTPUT_SEISMOS   = 5000:" DATA/Par_file
sedi "s:^NTSTEP_BETWEEN_OUTPUT_IMAGES .*:NTSTEP_BETWEEN_OUTPUT_IMAGES    = 5000:" DATA/Par_file
sedi "s:^output_color_image .*:output_color_image              = .false.:" DATA/Par_file
sedi "s:^output_postscript_snapshot .*:output_postscript_snapshot      = .false.:" DATA/Par_file
sedi "s:^output_grid_Gnuplot .*:output_grid_Gnuplot             = .false.:" DATA/Par_file
# the vertical element count lives in the interfaces file
sedi "s:^200$:100:" DATA/interfaces_elastic.dat

if [[ $? -ne 0 ]]; then
  echo "setup failed, please check..." >> $testdir/results.log
  exit 1
fi

###################################################
echo "step 1: mesher, picking the coupling elements" >> $testdir/results.log
###################################################
sedi "s:^COUPLING_IN .*:COUPLING_IN                     = .true.:" DATA/Par_file
run_solver xmeshfem2D
if [[ $? -ne 0 ]]; then echo "mesher failed, please check..." >> $testdir/results.log; exit 1; fi

if [ ! -e OUTPUT_FILES/externalsource.txt ]; then
  echo "no externalsource.txt written, please check..." >> $testdir/results.log
  exit 1
fi
cp OUTPUT_FILES/externalsource.txt ./externalsource.txt

# global element id of each coupling element, in the order the file lists them
awk -F',' '!/^ *#/ && NF >= 4 { gsub(/ /,"",$1); print $1 }' \
    ./externalsource.txt > ./ext_iele.txt
NCE=`wc -l < ./ext_iele.txt`
echo "  coupling elements: $NCE" >> $testdir/results.log
if [ "$NCE" -lt 10 ]; then
  echo "too few coupling elements ($NCE), please check dR_ext..." >> $testdir/results.log
  exit 1
fi

###################################################
echo "step 2: receivers on the coupling elements, plus the validation ring" >> $testdir/results.log
###################################################
# receivers 0 .. NCE-1 sit on the coupling elements, in externalsource.txt order;
# the validation receivers follow, outside the coupling surface
awk -F',' '!/^ *#/ && NF >= 4 { printf "S%04d  AA, %20.8f, %20.8f 0.0 0.0\n", n++, $3, $4 }' \
    ./externalsource.txt > DATA/STATIONS
awk -v n="$NCE" -v r="$R_VALIDATE" -v azs="$AZIMUTHS" 'BEGIN{
  m = split(azs, a, " ");
  pi = atan2(0,-1);
  for (i = 1; i <= m; i++) {
    # azimuth measured clockwise from north, north = +z, east = +x
    printf "S%04d  AA, %20.8f, %20.8f 0.0 0.0\n", n+i-1, r*sin(a[i]*pi/180.0), r*cos(a[i]*pi/180.0);
  }
}' >> DATA/STATIONS
NVAL=`echo $AZIMUTHS | wc -w`
echo "  receivers: $NCE on the coupling surface + $NVAL for validation" >> $testdir/results.log

###################################################
echo "step 3: point-source run (the reference)" >> $testdir/results.log
###################################################
sedi "s:^COUPLING_IN .*:COUPLING_IN                     = .false.:" DATA/Par_file
sedi "s:^factor .*:factor                          = 1.0d10:" DATA/SOURCE
rm -rf OUTPUT_FILES; mkdir -p OUTPUT_FILES
run_solver xmeshfem2D
if [[ $? -ne 0 ]]; then echo "mesher failed, please check..." >> $testdir/results.log; exit 1; fi
run_solver xspecfem2D
if [[ $? -ne 0 ]]; then echo "point-source run failed, please check..." >> $testdir/results.log; exit 1; fi
mv OUTPUT_FILES OUTPUT_FILES_P

###################################################
echo "step 4: turn those seismograms into the injected source" >> $testdir/results.log
###################################################
mkdir -p extsource
i=0
while read -r iele; do
  fx=`printf "OUTPUT_FILES_P/AA.S%04d.BXX.sema" $i`
  fz=`printf "OUTPUT_FILES_P/AA.S%04d.BXZ.sema" $i`
  if [ ! -e "$fx" ] || [ ! -e "$fz" ]; then
    echo "missing seismogram $fx or $fz, please check..." >> $testdir/results.log
    exit 1
  fi
  out=`printf "extsource/EXT%08d.dat" $iele`
  paste "$fx" "$fz" | awk '{ printf "%20.8f, %20.8e, %20.8e\n", $1, $2, $4 }' > "$out"
  i=`expr $i + 1`
done < ./ext_iele.txt
echo "  wrote $i EXT files" >> $testdir/results.log

###################################################
echo "step 5: coupling run, built-in source switched off" >> $testdir/results.log
###################################################
sedi "s:^COUPLING_IN .*:COUPLING_IN                     = .true.:" DATA/Par_file
sedi "s:^factor .*:factor                          = 0.0:" DATA/SOURCE
mkdir -p OUTPUT_FILES
run_solver xmeshfem2D
if [[ $? -ne 0 ]]; then echo "mesher failed, please check..." >> $testdir/results.log; exit 1; fi
run_solver xspecfem2D
if [[ $? -ne 0 ]]; then echo "coupling run failed, please check..." >> $testdir/results.log; exit 1; fi
mv OUTPUT_FILES OUTPUT_FILES_C

###################################################
echo "step 6: the two must agree outside the coupling surface" >> $testdir/results.log
###################################################
status=0
i=$NCE
for az in $AZIMUTHS; do
  for comp in BXX BXZ; do
    p=`printf "OUTPUT_FILES_P/AA.S%04d.%s.sema" $i $comp`
    c=`printf "OUTPUT_FILES_C/AA.S%04d.%s.sema" $i $comp`
    if [ ! -e "$p" ] || [ ! -e "$c" ]; then
      echo "  missing $p or $c" >> $testdir/results.log
      status=1
      continue
    fi
    r=`paste "$p" "$c" | awk '
      { n++; x=$2; y=$4; sx+=x; sy+=y; sxx+=x*x; syy+=y*y; sxy+=x*y }
      END {
        if (n == 0) { print "0.0"; exit }
        cov = sxy/n - (sx/n)*(sy/n);
        vx  = sxx/n - (sx/n)*(sx/n);
        vy  = syy/n - (sy/n)*(sy/n);
        if (vx <= 0 || vy <= 0) { print "0.0"; exit }
        printf "%.5f", cov/sqrt(vx*vy);
      }'`
    ok=`awk -v r="$r" -v t="$MIN_CORRELATION" 'BEGIN{ print (r >= t) ? 1 : 0 }'`
    if [ "$ok" -eq 1 ]; then
      echo "  azimuth ${az}deg $comp : correlation $r  good" >> $testdir/results.log
    else
      echo "  azimuth ${az}deg $comp : correlation $r  FAILED (< $MIN_CORRELATION)" >> $testdir/results.log
      status=1
    fi
  done
  i=`expr $i + 1`
done

if [ $status -ne 0 ]; then
  echo >> $testdir/results.log
  echo "coupling does not reproduce the point source, please check..." >> $testdir/results.log
  exit 1
fi

# cleanup
rm -rf ./DATA ./OUTPUT_FILES ./OUTPUT_FILES_P ./OUTPUT_FILES_C ./extsource
rm -f ./externalsource.txt ./ext_iele.txt

echo >> $testdir/results.log
echo "successful coupling validation" >> $testdir/results.log
