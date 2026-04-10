# References:
# - Tutorial 34: Solvation with 3D-RISM
#   https://ambermd.org/tutorials/advanced/tutorial34/index.html
# - Tutorial 40: 1D-RISM and 3D-RISM
#   https://ambermd.org/tutorials/advanced/tutorial40/index.php

import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from ..logger import generate_logger

LOGGER = generate_logger(__name__)

# sander 用 3D-RISM 入力テンプレート
SANDER_RISM_INPUT_TEMPLATE = """\
&cntrl
  ntx=1, nstlim=0, irism=1,
/
&rism
  closure='{closure}',
  grdspc={grdspc},{grdspc},{grdspc},
  tolerance={tolerance},
  buffer={buffer},
  solvcut={solvcut},
  mdiis_del=0.7,
  mdiis_nvec=5,
  maxstep=10000,
  npropagate=5,
  verbose=2,
  apply_rism_force=0,
  volfmt='dx',
  ntwrism=1,
/
"""


def add_subcmd(subparsers):
    """
    mdtbx place_solvent -p leap.parm7 -x leap.rst7 -o solvent_placed.pdb
    """
    parser = subparsers.add_parser(
        "place_solvent",
        help="Place solvent molecules using 3D-RISM",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-p",
        "--prmtop",
        required=True,
        type=str,
        help="AMBER topology file (.parm7)",
    )

    parser.add_argument(
        "-x",
        "--coord",
        required=True,
        type=str,
        help="AMBER coordinate file (.rst7)",
    )

    parser.add_argument(
        "-o",
        "--output",
        default="solvent_placed.pdb",
        type=str,
        help="Output PDB file",
    )

    parser.add_argument(
        "--xvv",
        default=None,
        type=str,
        help="Pre-computed solvent susceptibility file (.xvv). "
        "If not provided, rism1d is run to generate one.",
    )

    parser.add_argument(
        "--solvent",
        default="water",
        choices=["water"],
        help="Solvent type",
    )

    parser.add_argument(
        "--solvent-model",
        default="SPC",
        type=str,
        help="Solvent model for 1D-RISM (.mdl file stem in $AMBERHOME/dat/rism1d/mdl/)",
    )

    parser.add_argument(
        "--temperature",
        default=300.0,
        type=float,
        help="Temperature [K]",
    )

    parser.add_argument(
        "--closure",
        default="kh",
        choices=["kh", "hnc", "pse2", "pse3"],
        help="Closure approximation",
    )

    parser.add_argument(
        "--grdspc",
        default=0.5,
        type=float,
        help="Grid spacing [Å]",
    )

    parser.add_argument(
        "--tolerance",
        default=1e-5,
        type=float,
        help="Convergence tolerance for 3D-RISM",
    )

    parser.add_argument(
        "--buffer",
        default=14.0,
        type=float,
        help="Buffer distance around solute [Å]",
    )

    parser.add_argument(
        "--solvcut",
        default=14.0,
        type=float,
        help="Solvent cutoff distance [Å]",
    )

    parser.add_argument(
        "--threshold",
        default=1.5,
        type=float,
        help="Density threshold for peak extraction (g(r) value)",
    )

    parser.add_argument(
        "--exclusion-radius",
        default=2.6,
        type=float,
        help="Minimum distance between placed solvent sites [Å]. "
        "Approximate diameter of a water molecule.",
    )

    parser.add_argument(
        "--max-sites",
        default=None,
        type=int,
        help="Maximum number of solvent sites to place. "
        "If not set, all sites above threshold are placed.",
    )

    parser.add_argument(
        "--use-sander",
        action="store_true",
        help="Use sander interface instead of rism3d.snglpnt",
    )

    parser.add_argument(
        "--keepfiles",
        action="store_true",
        help="Keep intermediate files",
    )

    parser.set_defaults(func=run)


# ---------------------------------------------------------------------------
#  1D-RISM: xvv ファイル生成
# ---------------------------------------------------------------------------


def _run_rism1d(solvent_model, temperature, workdir):
    """1D-RISM を実行して xvv ファイルを生成する。

    $AMBERHOME/dat/rism1d/mdl/<solvent_model>.mdl が必要。
    Returns:
        xvv_path: 生成された xvv ファイルのパス
    """
    amberhome = os.environ.get("AMBERHOME", "")
    mdl_path = os.path.join(amberhome, "dat", "rism1d", "mdl", f"{solvent_model}.mdl")
    if not os.path.exists(mdl_path):
        raise FileNotFoundError(
            f"Solvent model file not found: {mdl_path}. "
            f"Check $AMBERHOME and --solvent-model."
        )

    xvv_stem = f"{solvent_model}_{temperature:.2f}"
    xvv_path = os.path.join(workdir, f"{xvv_stem}.xvv")

    # rism1d 入力ファイル
    inp_content = (
        f"&PARAMETERS\n"
        f"  THEORY='DRISM', CLOSURE='KH',\n"
        f"  NR=16384, DR=0.025,\n"
        f"  OUTLST='xvv',\n"
        f"  NIS=20, DESSION=0.5, MDIIS_NVEC=20, MDIIS_DEL=0.3,\n"
        f"  TOLERANCE=1.0e-12,\n"
        f"  SMEAR=1, APTS=0.2,\n"
        f"  TEMPER={temperature},\n"
        f"/\n"
        f"  {solvent_model}\n"
    )
    inp_path = os.path.join(workdir, f"{xvv_stem}.inp")
    with open(inp_path, "w") as f:
        f.write(inp_content)

    rism1d_cmd = f"rism1d {xvv_stem} > {xvv_stem}.out 2>&1"
    LOGGER.info(f"Running 1D-RISM to generate xvv file ({solvent_model}) ...")
    subprocess.run(rism1d_cmd, shell=True, check=True, cwd=workdir)

    if not os.path.exists(xvv_path):
        raise RuntimeError(
            f"1D-RISM did not produce {xvv_path}. Check the output for errors."
        )

    LOGGER.info(f"Generated xvv file: {xvv_path}")
    return xvv_path


# ---------------------------------------------------------------------------
#  3D-RISM 実行
# ---------------------------------------------------------------------------


def _run_rism3d_snglpnt(prmtop, coord, xvv_path, args, workdir):
    """rism3d.snglpnt コマンドラインインターフェースで 3D-RISM を実行する。"""

    prmtop_abs = os.path.abspath(prmtop)
    coord_abs = os.path.abspath(coord)
    xvv_abs = os.path.abspath(xvv_path)
    prmtop_stem = Path(prmtop).stem

    guv_prefix = os.path.join(workdir, prmtop_stem)

    rism_cmd = [
        "rism3d.snglpnt",
        "--prmtop",
        prmtop_abs,
        "--rst",
        coord_abs,
        "--xvv",
        xvv_abs,
        "--closure",
        args.closure,
        "--grdspc",
        f"{args.grdspc},{args.grdspc},{args.grdspc}",
        "--buffer",
        str(args.buffer),
        "--solvcut",
        str(args.solvcut),
        "--tolerance",
        str(args.tolerance),
        "--verbose",
        "2",
        "--ntwrism",
        "1",
        "--guv",
        guv_prefix,
    ]

    LOGGER.info("Running rism3d.snglpnt ...")
    LOGGER.info(f"  Command: {' '.join(rism_cmd)}")
    subprocess.run(rism_cmd, check=True, cwd=workdir)


def _run_sander_rism(prmtop, coord, args, workdir):
    """sander インターフェースで 3D-RISM を実行する。"""

    prmtop_abs = os.path.abspath(prmtop)
    coord_abs = os.path.abspath(coord)

    mdin_content = SANDER_RISM_INPUT_TEMPLATE.format(
        closure=args.closure,
        grdspc=args.grdspc,
        tolerance=args.tolerance,
        buffer=args.buffer,
        solvcut=args.solvcut,
    )
    mdin_path = os.path.join(workdir, "mdin.rism")
    with open(mdin_path, "w") as f:
        f.write(mdin_content)

    sander_cmd = [
        "sander",
        "-O",
        "-i",
        mdin_path,
        "-o",
        os.path.join(workdir, "mdout"),
        "-p",
        prmtop_abs,
        "-c",
        coord_abs,
    ]

    LOGGER.info("Running sander with 3D-RISM ...")
    subprocess.run(sander_cmd, check=True, cwd=workdir)


# ---------------------------------------------------------------------------
#  DX ファイル読み込み
# ---------------------------------------------------------------------------


def _parse_dx(dx_path):
    """OpenDX 形式の密度グリッドを ndarray として読み込む。

    Returns:
        data: (nx, ny, nz) ndarray — g(r) 分布関数値
        origin: (3,) ndarray — グリッド原点座標 [Å]
        delta: (3, 3) ndarray — 各軸のグリッド刻み幅ベクトル (行ごと)
    """
    origin = None
    delta = []
    counts = None
    data_values = []

    with open(dx_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            if line.startswith("object 1"):
                parts = line.split()
                counts = tuple(int(x) for x in parts[-3:])
            elif line.startswith("origin"):
                parts = line.split()
                origin = np.array([float(x) for x in parts[1:4]])
            elif line.startswith("delta"):
                parts = line.split()
                delta.append([float(x) for x in parts[1:4]])
            elif line.startswith("object") or line.startswith("attribute"):
                continue
            elif line.startswith("component"):
                continue
            else:
                # データ行
                data_values.extend(float(v) for v in line.split())

    if counts is None:
        raise ValueError(f"Could not parse grid dimensions from {dx_path}")
    if origin is None:
        raise ValueError(f"Could not parse origin from {dx_path}")
    if len(delta) != 3:
        raise ValueError(f"Expected 3 delta vectors, got {len(delta)}")

    data = np.array(data_values).reshape(counts)
    delta = np.array(delta)

    return data, origin, delta


def _grid_to_cartesian(indices, origin, delta):
    """グリッドインデックス (N, 3) を直交座標 (N, 3) に変換する。

    OpenDX の delta 行列は行ごとに各軸方向の刻み幅ベクトルを持つ。
    座標 = origin + i * delta[0] + j * delta[1] + k * delta[2]
         = origin + indices @ delta
    """
    return origin + indices @ delta


# ---------------------------------------------------------------------------
#  Placevent 風グリーディピーク抽出
# ---------------------------------------------------------------------------


def _extract_peaks_greedy(
    data, origin, delta, threshold, exclusion_radius, max_sites=None
):
    """Placevent アルゴリズムに基づくグリーディな溶媒サイト抽出。

    1. g(r) > threshold のグリッド点をすべて候補とする
    2. g(r) 値の降順にソートする
    3. 最も高い g(r) の点を溶媒サイトとして採用し、
       exclusion_radius 以内の他の候補をすべて除外する
    4. 残りの候補で最も高い g(r) の点を次のサイトとする
    5. max_sites に達するか候補がなくなるまで繰り返す

    Returns:
        coords: (M, 3) ndarray — 溶媒サイト座標 [Å]
        gvalues: (M,) ndarray — 各サイトの g(r) 値
    """
    # 閾値を超えるグリッド点のインデックスと値を取得
    indices = np.argwhere(data > threshold)
    if len(indices) == 0:
        LOGGER.warning(
            f"No grid points exceed threshold g(r) > {threshold}. "
            "Try lowering --threshold."
        )
        return np.empty((0, 3)), np.empty(0)

    values = data[indices[:, 0], indices[:, 1], indices[:, 2]]

    # g(r) 降順でソート
    order = np.argsort(-values)
    indices = indices[order]
    values = values[order]

    # 座標に変換
    all_coords = _grid_to_cartesian(indices.astype(float), origin, delta)

    # グリーディ選択
    placed_coords = []
    placed_gvalues = []
    used = np.zeros(len(all_coords), dtype=bool)
    excl_sq = exclusion_radius**2

    for i in range(len(all_coords)):
        if used[i]:
            continue

        coord_i = all_coords[i]
        placed_coords.append(coord_i)
        placed_gvalues.append(values[i])

        if max_sites is not None and len(placed_coords) >= max_sites:
            break

        # この点から exclusion_radius 以内の候補を除外
        remaining = np.where(~used)[0]
        remaining = remaining[remaining > i]
        if len(remaining) > 0:
            diff = all_coords[remaining] - coord_i
            dist_sq = np.sum(diff**2, axis=1)
            too_close = remaining[dist_sq < excl_sq]
            used[too_close] = True

    coords = np.array(placed_coords)
    gvalues = np.array(placed_gvalues)

    LOGGER.info(
        f"Extracted {len(coords)} solvent sites "
        f"(threshold={threshold}, exclusion_radius={exclusion_radius} Å)"
    )
    return coords, gvalues


# ---------------------------------------------------------------------------
#  PDB 出力
# ---------------------------------------------------------------------------


def _write_pdb(coords, gvalues, solvent, output_path):
    """溶媒サイト座標を PDB 形式で書き出す。

    occupancy に g(r) の初期値を、B-factor に配置順の g(r) を記録する。
    原子番号・残基番号が PDB フォーマットの上限を超える場合は
    モジュロで折り返す。
    """
    atom_name = "O" if solvent == "water" else "X"
    resname = "WAT" if solvent == "water" else "SOL"

    with open(output_path, "w") as f:
        f.write(
            "REMARK   Generated by place_solvent (3D-RISM)\n"
            f"REMARK   {len(coords)} solvent sites placed\n"
        )
        for i, (coord, gval) in enumerate(zip(coords, gvalues), start=1):
            x, y, z = coord
            serial = i % 100000  # ATOM serial は 5 桁まで
            resseq = i % 10000  # 残基番号は 4 桁まで
            f.write(
                f"HETATM{serial:5d}  {atom_name:<3s} {resname} A"
                f"{resseq:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}"
                f"{1.00:6.2f}{gval:6.2f}\n"
            )
        f.write("END\n")

    LOGGER.info(f"{len(coords)} solvent sites written to {output_path}")


# ---------------------------------------------------------------------------
#  DX ファイル検索
# ---------------------------------------------------------------------------


def _find_oxygen_dx(workdir, prmtop_stem, closure):
    """3D-RISM が出力した酸素密度の .dx ファイルを見つける。

    rism3d.snglpnt --guv prefix の場合:
        prefix.O.1.dx (guv 出力)
    sander の場合:
        <stem>.<closure>.O.0.dx (guv 出力, 0-indexed)

    いずれも酸素サイト "O" を含むファイルを探す。
    """
    workpath = Path(workdir)

    # guv 出力パターン（rism3d.snglpnt）
    # prefix.O.1.dx が典型的
    candidates = sorted(workpath.glob(f"{prmtop_stem}*O*.dx"))

    if not candidates:
        # sander 出力パターン
        candidates = sorted(workpath.glob(f"*{closure}*O*.dx"))

    if not candidates:
        # フォールバック: 任意の酸素を含む dx
        candidates = sorted(workpath.glob("*O*.dx"))

    if not candidates:
        # 最終手段: 全 dx ファイル
        candidates = sorted(workpath.glob("*.dx"))

    if not candidates:
        raise FileNotFoundError(
            f"No .dx output files found in {workdir}. "
            "3D-RISM calculation may have failed."
        )

    # 水の酸素に最も適合するファイルを選択
    # "guv" や "g" を含むもの (分布関数) を優先し、
    # "cuv" (直接相関) や "huv" (間接相関) を避ける
    for cand in candidates:
        name = cand.name.lower()
        if "cuv" in name or "huv" in name or "uuv" in name:
            continue
        return cand

    # すべて除外された場合は最初のものを使う
    return candidates[0]


# ---------------------------------------------------------------------------
#  メイン
# ---------------------------------------------------------------------------


def run(args):
    prmtop_stem = Path(args.prmtop).stem

    # 作業ディレクトリ
    workdir = tempfile.mkdtemp(prefix="rism3d_")
    LOGGER.info(f"Working directory: {workdir}")

    try:
        # ---- 1. xvv ファイル準備 ----
        if args.xvv is not None:
            xvv_path = os.path.abspath(args.xvv)
            if not os.path.exists(xvv_path):
                LOGGER.error(f"xvv file not found: {xvv_path}")
                return
            LOGGER.info(f"Using provided xvv file: {xvv_path}")
        else:
            xvv_path = _run_rism1d(args.solvent_model, args.temperature, workdir)

        # ---- 2. 3D-RISM 実行 ----
        if args.use_sander:
            _run_sander_rism(args.prmtop, args.coord, args, workdir)
        else:
            _run_rism3d_snglpnt(args.prmtop, args.coord, xvv_path, args, workdir)

        # ---- 3. 酸素密度 DX ファイル読み込み ----
        dx_path = _find_oxygen_dx(workdir, prmtop_stem, args.closure)
        LOGGER.info(f"Reading density from {dx_path}")

        data, origin, delta = _parse_dx(str(dx_path))
        LOGGER.info(
            f"Grid dimensions: {data.shape}, "
            f"origin: {origin}, "
            f"spacing: {np.diag(delta)}"
        )

        # ---- 4. ピーク抽出 (Placevent 風グリーディ) ----
        coords, gvalues = _extract_peaks_greedy(
            data,
            origin,
            delta,
            threshold=args.threshold,
            exclusion_radius=args.exclusion_radius,
            max_sites=args.max_sites,
        )

        if len(coords) == 0:
            LOGGER.error("No solvent sites found.")
            return

        # ---- 5. PDB 出力 ----
        output_path = os.path.abspath(args.output)
        _write_pdb(coords, gvalues, args.solvent, output_path)

    finally:
        if not args.keepfiles:
            shutil.rmtree(workdir, ignore_errors=True)
            LOGGER.info("Intermediate files removed")
        else:
            LOGGER.info(f"Intermediate files kept in {workdir}")
