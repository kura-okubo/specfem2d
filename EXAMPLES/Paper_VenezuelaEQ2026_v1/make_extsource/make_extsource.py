#!/usr/bin/env python3
"""HOSS のカップリングセンサー出力から SPECFEM2D の EXT ファイルを作る。

    python3 make_extsource.py sensors_v17_start.tar.gz
    python3 make_extsource.py /path/to/extracted/sensors/

センサーは速度を出す。SPECFEM2D の外部震源は加速度を要求する（compute_ext_source.F90
の accel_elastic = extsource）。素直に差分を取ると高周波が増幅されるので、
微分・低域通過・内挿を 1 回の FFT でまとめて行う。

  1. 端をコサインテーパで 0 に落とす。HOSS は 100 s で終わるが速度はまだ残って
     いるので、切りっぱなしだと加速度に段差パルスが出る。
  2. スペクトルで i*omega を掛けて微分し、同時に FC で低域通過。メッシュが
     伝播できない帯域（solver のログの maximum dominant source frequency）は
     先に落とす。
  3. 周波数領域のゼロ詰めで HOSS の 20 Hz から SPECFEM2D の DT へ帯域制限内挿。

要素とセンサーの対応は venezuela_coupling_sensors.csv（ノートブックが書く
カップリングマニフェスト）から取る。Sensor.input の行番号は on-fault の本数が
変わるとずれるので、ここを固定値で書いてはいけない。
"""
import argparse, os, sys, tarfile, tempfile, time
import numpy as np
import pandas as pd

FC, FSTOP = 0.5, 0.75                 # Hz, 低域通過（平坦端 / ゼロ端）
TAPER_HEAD, TAPER_TAIL = 1.0, 15.0    # s
HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.dirname(HERE)
DEFAULT_MANIFEST = os.path.expanduser(
    "~/Library/CloudStorage/Dropbox/NIED_RESEARCH/VenezuelaEQ_Paper_v1/"
    "StressState/data/venezuela_coupling_sensors.csv")


def taper(n, dt, head, tail):
    w = np.ones(n)
    nh, nt = int(head / dt), int(tail / dt)
    w[:nh] = .5 * (1 - np.cos(np.pi * np.arange(nh) / nh))
    w[n - nt:] = .5 * (1 + np.cos(np.pi * np.arange(nt) / nt))
    return w


def velocity_to_accel(v, dt_in, dt_out, n_out):
    nt = v.shape[1]
    w = taper(nt, dt_in, TAPER_HEAD, TAPER_TAIL)
    n0 = 1 << int(np.ceil(np.log2(nt * 1.25)))       # ゼロ詰めで巻き込みを防ぐ
    vp = np.zeros((v.shape[0], n0)); vp[:, :nt] = v * w
    up = int(round(dt_in / dt_out))
    if abs(dt_in / dt_out - up) > 1e-9:
        sys.exit(f"dt_in/dt_out = {dt_in/dt_out} is not an integer")
    f = np.fft.rfftfreq(n0, dt_in)
    H = np.ones_like(f)
    band = (f > FC) & (f < FSTOP)
    H[band] = .5 * (1 + np.cos(np.pi * (f[band] - FC) / (FSTOP - FC)))
    H[f >= FSTOP] = 0.0
    A = np.fft.rfft(vp, axis=1) * (2j * np.pi * f) * H
    B = np.zeros((v.shape[0], n0 * up // 2 + 1), dtype=complex)
    B[:, :A.shape[1]] = A
    return (np.fft.irfft(B, n=n0 * up, axis=1) * up)[:, :n_out]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sensors", help="sensors tar.gz, or a directory of sensorNNNN files")
    ap.add_argument("--manifest", default=DEFAULT_MANIFEST)
    ap.add_argument("--case", default=CASE)
    ap.add_argument("--seconds", type=float, default=105.0,
                    help="how much of the trace to write; the solver zero-fills the rest")
    a = ap.parse_args()

    par = open(os.path.join(a.case, "DATA", "Par_file")).read()
    import re
    dt_out = float(re.search(r"^\s*DT\s*=\s*([0-9.eEdD+-]+)", par, re.M)
                   .group(1).lower().replace("d", "e"))
    nstep = int(re.search(r"^\s*NSTEP\s*=\s*(\d+)", par, re.M).group(1))
    n_out = min(int(round(a.seconds / dt_out)), nstep)
    print(f"Par_file: DT = {dt_out} s, NSTEP = {nstep} -> writing {n_out} samples "
          f"({n_out*dt_out:.0f} s); the solver zero-fills to {nstep*dt_out:.0f} s")

    man = pd.read_csv(a.manifest)
    print(f"manifest: {len(man)} coupling elements, "
          f"sensor_index {man.sensor_index.min()}..{man.sensor_index.max()}")

    tmp = None
    src = a.sensors
    if os.path.isfile(src):
        tmp = tempfile.mkdtemp()
        print(f"extracting {src} ...")
        with tarfile.open(src) as tf:
            tf.extractall(tmp)
        src = tmp

    # the archive is usually made from the stage's parent, so the files sit in
    # start/ rather than at the top; walk down while there is exactly one
    # directory and no sensor files here
    while not any(f.startswith("sensor") for f in os.listdir(src)):
        sub = [f for f in os.listdir(src) if os.path.isdir(os.path.join(src, f))]
        if len(sub) != 1:
            sys.exit(f"cannot find the sensor files under {src}")
        src = os.path.join(src, sub[0])
        print(f"  descending into {os.path.basename(src)}/")

    # Sensor.input の行番号は 1 始まり、HOSS のファイル名は 0 始まり
    cols = ["t", "cx", "cy", "vx", "vy", "Cxx", "Cyy", "Cxy"]
    VX = VY = None
    t0 = time.time()
    for k, idx in enumerate(man.sensor_index.to_numpy()):
        p = os.path.join(src, "sensor%04d" % (idx - 1))
        d = pd.read_csv(p, comment="!", header=None, names=cols).to_numpy()
        if VX is None:
            t = d[:, 0]
            VX = np.empty((len(man), len(t))); VY = np.empty_like(VX)
        VX[k], VY[k] = d[:, 3], d[:, 4]
        if k % 200 == 0:
            print(f"  read {k}/{len(man)} ({time.time()-t0:.0f} s)")
    dt_in = t[1] - t[0]
    print(f"sensors: {len(man)} x {len(t)} samples, dt = {dt_in} s "
          f"({1/dt_in:.0f} Hz), t = {t[0]:.1f}..{t[-1]:.1f} s")
    print(f"  |v| max {np.hypot(VX, VY).max():.4f} m/s, "
          f"|v| at the last sample {np.hypot(VX[:, -1], VY[:, -1]).max():.3e} m/s "
          f"(tapered to zero over the last {TAPER_TAIL:.0f} s)")

    ax = velocity_to_accel(VX, dt_in, dt_out, n_out)
    az = velocity_to_accel(VY, dt_in, dt_out, n_out)   # HOSS の y(北) = SPECFEM2D の z
    if not (np.isfinite(ax).all() and np.isfinite(az).all()):
        sys.exit("non-finite acceleration; the solver would abort on the NaN check")
    print(f"acceleration: |a| max {np.hypot(ax, az).max():.4e} m/s^2")

    # 積分して速度に戻るか（処理が壊れていないことの確認）
    kmax = int(np.argmax(np.hypot(VX, VY).max(axis=1)))
    vint = np.cumsum(ax[kmax]) * dt_out
    vref = np.interp(np.arange(n_out) * dt_out, t,
                     VX[kmax] * taper(len(t), dt_in, TAPER_HEAD, TAPER_TAIL))
    m = np.arange(n_out) * dt_out <= t[-1]
    print(f"  check: integrating back gives r = {np.corrcoef(vint[m], vref[m])[0,1]:.6f}, "
          f"amplitude ratio {np.abs(vint[m]).max()/np.abs(vref[m]).max():.4f}")

    out = os.path.join(a.case, "extsource")
    os.makedirs(out, exist_ok=True)
    tcol = np.char.mod("%.5f", np.arange(n_out) * dt_out)
    t0 = time.time()
    for k, e in enumerate(man.iele.to_numpy()):
        body = "\n".join("%s %.6e %.6e" % r for r in zip(tcol, ax[k], az[k]))
        open(os.path.join(out, "EXT%08d.dat" % e), "w").write(body + "\n")
    print(f"wrote {len(man)} files to {out} ({time.time()-t0:.0f} s)")

    # 実機のメッシャが出した要素リストと突き合わせる（これが正）
    ext = os.path.join(a.case, "OUTPUT_FILES_grid", "externalsource.txt")
    if os.path.exists(ext):
        ids = [int(l.split(",")[0]) for l in open(ext)
               if l.strip() and not l.strip().startswith("#")]
        miss = set(ids) - set(man.iele)
        extra = set(man.iele) - set(ids)
        print(f"externalsource.txt: {len(ids)} elements, "
              f"{len(miss)} missing here, {len(extra)} extra")
        if miss or extra:
            sys.exit("element list disagrees with the mesher -- do not run")
    else:
        print(f"NOTE: {ext} not found; run run_xmeshfem2D.sh to cross-check")
    if tmp:
        import shutil; shutil.rmtree(tmp)


if __name__ == "__main__":
    main()
