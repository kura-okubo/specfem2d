#!/bin/bash
# point-source trial run (step 5); coupling off
set -u
echo "running: `date`"

# scalar moment per unit out-of-plane length, M0 / W
FACTOR=5.127799e+15
sed -i.bak "s/^COUPLING_IN .*/COUPLING_IN                     = .false./" DATA/Par_file
sed -i.bak "s/^factor .*/factor                          = $FACTOR/" DATA/SOURCE
rm -f DATA/*.bak

rm -rf OUTPUT_FILES_P
mkdir -p OUTPUT_FILES
rm -rf OUTPUT_FILES/*

rm -f xmeshfem2D xspecfem2D
ln -s ../../bin/xmeshfem2D ./
ln -s ../../bin/xspecfem2D ./

cp DATA/Par_file DATA/SOURCE DATA/STATIONS OUTPUT_FILES/

NPROC=`grep ^NPROC DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`

echo
echo "running mesher on $NPROC processors..."
echo
mpirun -np $NPROC ./xmeshfem2D
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "running solver on $NPROC processors..."
echo
mpirun -np $NPROC ./xspecfem2D
if [[ $? -ne 0 ]]; then exit 1; fi

# --- mesh check, as PNG -----------------------------------------------------
# xmeshfem2D only offers PostScript, which macOS cannot preview. Render it.
gnuplot OUTPUT_FILES/plot_gridfile.gnu < /dev/null || true

if [ -f OUTPUT_FILES/gridfile.ps ]; then
  # Faithful conversion of what the mesher drew. The PostScript is landscape
  # on a portrait page and gs does not honour %%Orientation, so rotate after.
  gs -dSAFER -dBATCH -dNOPAUSE -dQUIET -sDEVICE=png16m -r200      -dTextAlphaBits=4 -dGraphicsAlphaBits=4      -sOutputFile=OUTPUT_FILES/gridfile.png OUTPUT_FILES/gridfile.ps      && sips -r 90 OUTPUT_FILES/gridfile.png >/dev/null 2>&1
fi

# That one goes solid black as soon as the mesh is fine (at h = 1.5 km the
# element outlines fill the page), so also draw readable versions from the
# same data: the mesh in pale grey, the coupling elements in red.
if [ -f OUTPUT_FILES/gridfile.gnu ]; then
  # zoomed on the coupling band -- the one to actually check
  gnuplot <<'GNUPLOT' || true
set term pngcairo size 3200,1308 font "Helvetica,18"
set output "OUTPUT_FILES/gridfile_coupling.png"
set size ratio -1
set title "Coupling elements (red) on the SPECFEM2D mesh (grey), h = 1.5 km"
set xlabel "x [m]  (UTM 19N easting - origin)"
set ylabel "z [m]  (UTM 19N northing - origin)"
set xrange [-100000.0:250000.0]
set yrange [-55500.0:68000.0]
set loadpath "./OUTPUT_FILES/"
load "gridfile_externalsource.gnu"
plot "gridfile.gnu" title "" w l lc rgb "#C8C8C8" lw 0.4
GNUPLOT

  # whole domain, for context
  gnuplot <<'GNUPLOT' || true
set term pngcairo size 3200,1900 font "Helvetica,18"
set output "OUTPUT_FILES/gridfile_coupling_full.png"
set size ratio -1
set title "Whole SPECFEM2D domain; coupling elements in red"
set xlabel "x [m]  (UTM 19N easting - origin)"
set ylabel "z [m]  (UTM 19N northing - origin)"
set loadpath "./OUTPUT_FILES/"
load "gridfile_externalsource.gnu"
plot "gridfile.gnu" title "" w l lc rgb "#E0E0E0" lw 0.2
GNUPLOT
fi

mv OUTPUT_FILES OUTPUT_FILES_P

echo
echo "see results in OUTPUT_FILES_P/"
echo `date`
