#!/usr/bin/env python3
"""SPECFEM2D のスナップショットに断層トレース等を重ねてパネル図にする。

    python3 plot_snapshots.py [OUTPUT_FILES_C] [out.png]

forward_image*.jpg は白枠付きで書かれるので、白でない領域を計算領域
(xmin..xmax, zmin..zmax) に対応づけて data 座標で重ね描きする。
"""
import glob, os, re, sys
import numpy as np
from PIL import Image
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
SRC  = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "OUTPUT_FILES_C")
OUT  = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "snapshots.png")

XMIN, XMAX = -600e3, 1062e3          # DATA/Par_file
ZMIN, ZMAX = -320e3, 571e3           # DATA/interfaces_venezuela.dat
DT   = 0.01
NPANEL = 12

def frame_extent(a):
    """白枠を除いた描画領域のピクセル範囲。"""
    nw = ~((a[:, :, 0] > 230) & (a[:, :, 1] > 230) & (a[:, :, 2] > 230))
    ys, xs = np.where(nw)
    return xs.min(), xs.max(), ys.min(), ys.max()

files = sorted(glob.glob(os.path.join(SRC, "forward_image*.jpg")))
if not files:
    sys.exit("no forward_image*.jpg in %s" % SRC)
steps = [int(re.search(r"(\d+)\.jpg", f).group(1)) for f in files]
pick  = np.linspace(0, len(files) - 1, min(NPANEL, len(files))).astype(int)

tr = {}
with open(os.path.join(HERE, "fault_traces.csv")) as fh:
    next(fh)
    for line in fh:
        n, x, z = line.split(",")
        tr.setdefault(n, []).append((float(x), float(z)))

# カップリング帯（Par_file の rec_*）
rec = dict(xmin=-75e3, zmin=-30.5e3, xmax=225e3, zmax=43e3)

sta = []
sf = os.path.join(HERE, "..", "DATA", "STATIONS")
if os.path.exists(sf):
    for line in open(sf):
        p = line.split()
        if len(p) >= 4: sta.append((float(p[2]), float(p[3])))
sta = np.array(sta)

ncol = 3; nrow = int(np.ceil(len(pick) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(5.6 * ncol, 3.4 * nrow),
                         sharex=True, sharey=True)
for ax, k in zip(np.atleast_1d(axes).ravel(), pick):
    a = np.asarray(Image.open(files[k]).convert("RGB"))
    x0, x1, y0, y1 = frame_extent(a)
    ax.imshow(a[y0:y1 + 1, x0:x1 + 1],
              extent=[XMIN / 1e3, XMAX / 1e3, ZMIN / 1e3, ZMAX / 1e3],
              origin="upper", interpolation="bilinear")
    for pts in tr.values():
        p = np.array(pts) / 1e3
        ax.plot(p[:, 0], p[:, 1], "-", c="lime", lw=1.0)
    ax.add_patch(plt.Rectangle((rec["xmin"] / 1e3, rec["zmin"] / 1e3),
                               (rec["xmax"] - rec["xmin"]) / 1e3,
                               (rec["zmax"] - rec["zmin"]) / 1e3,
                               fill=False, ec="cyan", lw=1.0, ls="--"))
    if len(sta): ax.plot(sta[:, 0] / 1e3, sta[:, 1] / 1e3, "w^", ms=3, mec="k", mew=.3)
    ax.set_title("t = %.0f s" % (steps[k] * DT), fontsize=10)
    ax.set_aspect("equal")
for ax in np.atleast_1d(axes).ravel()[len(pick):]: ax.axis("off")
for ax in np.atleast_2d(axes)[-1]: ax.set_xlabel("x [km]")
for ax in np.atleast_2d(axes)[:, 0]: ax.set_ylabel("z [km]")
fig.suptitle("Venezuela 2026 -- SPECFEM2D coupled run (HOSS v16 injected on the cyan box); "
             "velocity Vx", y=.995)
plt.tight_layout(rect=[0, 0, 1, .985])
fig.savefig(OUT, dpi=115, bbox_inches="tight")
print("wrote", OUT, "from", len(files), "frames in", SRC)
