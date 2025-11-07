#!/usr/bin/env python
# coding: utf-8

import pickle
import sys
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import List

import matplotlib as mpl
import numpy as np
import pandas as pd
import deeptime
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.spatial import ConvexHull
from scipy.interpolate import CubicSpline
import subprocess as sb


# input the following variables
n_trial_for_calc: List[int] = [i for i in range(1, 25 + 1)]
trial_root_directory: str = "./trials/"
feature_1d_directory: str = "./cvs_reshaped/comdist/"
feature_3d_directory: str = "./cvs_reshaped/comvec/"
output_directory: str = "./out_dist_msm"
show_picture: bool = False
T: float = 310
dt: int = 1
n_clusters_for_try_1d: List[int] = [
    20,
    50,
    100,
    500,
]
lags_for_try_1d: List[int] = [
    1,
    5,
    10,
    20,
    30,
    40,
    50,
]
n_clusters_1d_main = n_clusters_for_try_1d
lags_1d_main = lags_for_try_1d

n_clusters_for_try_3d: List[int] = [
    20,
    50,
    100,
    500,
    1000,
]
lags_for_try_3d: List[int] = [
    1,
    5,
    10,
    20,
    30,
    40,
    50,
]
n_clusters_3d_main = n_clusters_for_try_3d
lags_3d_main = lags_for_try_3d

cutoff: float = 12.5
nbins: int = 15
cmap: ListedColormap = mpl.colormaps.get_cmap("tab20")
do_volume_correction: bool = True
num_of_ligand: int = 1
box_size: float = 17 * 17 * 17
# definition for bound-state
lower_bound = 0  # 0 is recommended
upper_bound = 2

# definition for unbound-state
lower_unbound = 4
upper_unbound = cutoff  # params.cutoff is recommended


cmap_tab20 = mpl.colormaps.get_cmap("tab20")
cmap_tab20b = mpl.colormaps.get_cmap("tab20b")

colors_tab20_tab20b = np.vstack(
    (cmap_tab20(np.linspace(0, 1, 20)), cmap_tab20b(np.linspace(0, 1, 20)))
)
cmap_tab20_tab20b = ListedColormap(colors_tab20_tab20b)

cmap: ListedColormap = cmap_tab20_tab20b


rgb_list = [
    [1, 0, 0],  # 1: red
    [0.35, 0.35, 0.35],  # 2: gray
    [1, 0.5, 0],  # 3: orange
    [1, 1, 0],  # 4: yellow
    [0.5, 0.5, 0.2],  # 5: tan
    [0.6, 0.6, 0.6],  # 6: silver
    [0, 1, 0],  # 7: green
    [1, 1, 1],  # 8: white
    [1, 0.6, 0.6],  # 9: pink
    [0.25, 0.75, 0.75],  # 10: cyan
    [0.65, 0, 0.65],  # 11: purple
    [0.5, 0.9, 0.4],  # 12: lime
    [0.9, 0.4, 0.7],  # 13: mauve
    [0.5, 0.3, 0],  # 14: orche
    [0.5, 0.5, 0.75],  # 15: iceblue
    [0, 0, 0],  # 16: black
    [0.88, 0.97, 0.02],  # 17: yellow2
    [0.55, 0.9, 0.02],  # 18: yellow3
    [0, 0.9, 0.04],  # 19: green2
    [0, 0.9, 0.5],  # 20: green3
    [0, 0.88, 1],  # 21: cyan2
    [0, 0.76, 1],  # 22: cyan3
    [0.02, 0.38, 0.67],  # 23: blue2
    [0.01, 0.04, 0.93],  # 24: blue3
    [0.27, 0, 0.98],  # 25: violet
    [0.45, 0, 0.9],  # 26: violet2
    [0.9, 0, 0.9],  # 27: magenta
    [1, 0, 0.66],  # 28: magenta2
    [0.98, 0, 0.23],  # 29: red2
    [0.81, 0, 0],  # 30: red3
    [0.89, 0.35, 0],  # 31: orange2
    [0.96, 0.72, 0],  # 32: orange3
]

hex_colors = [mpl.colors.to_hex(rgb) for rgb in rgb_list]
vmd_cmap = ListedColormap(hex_colors)


# define the logger object and log file
logger = logging.getLogger(__name__)
log_file = Path(output_directory) / "logs/distance_deeptime.log"
log_file.parent.mkdir(parents=True, exist_ok=True)
logger.setLevel(logging.INFO)

# define the formatter
formatter = logging.Formatter(
    "%(asctime)s [%(levelname)s] (%(funcName)s) : %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

# file handler for logging (save to log file)
file_handler = logging.FileHandler(log_file)
file_handler.setLevel(logging.INFO)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

# stream handler for logging (print to console)
stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.INFO)
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)


# logging the error to both console and file
# a bit dirty way to handle unexpected errors, but it works
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        # KeyboardInterrupt is not an error
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    logger.error(
        "An unhandled exception occurred.",
        exc_info=(exc_type, exc_value, exc_traceback),
    )


sys.excepthook = handle_exception


@dataclass(frozen=True)
class Parameters:
    T: float
    trial_root_directory: str
    n_trial_for_calc: List[int]
    feature_1d_directory: str
    feature_3d_directory: str
    output_directory: str
    show_picture: bool
    dt: int
    n_clusters_for_try_1d: List[int]
    lags_for_try_1d: List[int]
    n_clusters_for_try_3d: List[int]
    lags_for_try_3d: List[int]
    cutoff: float
    nbins: int
    do_volume_correction: bool
    num_of_ligand: int
    box_size: float
    cmap: ListedColormap
    logger: logging.Logger = logger
    _RT: float = (
        -8.314 * T / 1000 / 4.184
    )  # constant for calculating free energy [kcal/mol/T]
    ligand_concentration: float = (
        num_of_ligand / box_size * (10 / 6.022)
    )  # [mol/L], ligand concentration used for koff calculation


# parameters
params = Parameters(
    T=T,
    trial_root_directory=trial_root_directory,
    n_trial_for_calc=n_trial_for_calc,
    feature_1d_directory=feature_1d_directory,
    feature_3d_directory=feature_3d_directory,
    output_directory=output_directory,
    show_picture=show_picture,
    dt=dt,
    n_clusters_for_try_1d=n_clusters_for_try_1d,
    lags_for_try_1d=lags_for_try_1d,
    n_clusters_for_try_3d=n_clusters_for_try_3d,
    lags_for_try_3d=lags_for_try_3d,
    cutoff=cutoff,
    nbins=nbins,
    do_volume_correction=do_volume_correction,
    num_of_ligand=num_of_ligand,
    box_size=box_size,
    cmap=cmap,  # use vmd_cmap if you want to use VMD color list
    logger=logger,
)

# create output directory
Path(params.output_directory).mkdir(parents=True, exist_ok=True)
Path(f"{params.output_directory}/images").mkdir(parents=True, exist_ok=True)
Path(f"{params.output_directory}/cluster_objs").mkdir(parents=True, exist_ok=True)
Path(f"{params.output_directory}/MSM_objs").mkdir(parents=True, exist_ok=True)
Path(f"{params.output_directory}/Count_objs").mkdir(parents=True, exist_ok=True)
Path(f"{params.output_directory}/result_csvs").mkdir(parents=True, exist_ok=True)
Path(f"{params.output_directory}/logs").mkdir(parents=True, exist_ok=True)


if not params.show_picture:
    mpl.use("Agg")


def plot_feature(params: Parameters, threshold: float) -> None:
    plt.figure(figsize=(20, 5))
    for trial in params.n_trial_for_calc:
        if not Path(f"{params.trial_root_directory}/trial{trial:03}/").exists():
            params.logger.info(f"trial{trial:03} was not found")
            continue
        else:
            res = sb.run(
                f"head -n 1 {params.trial_root_directory}/trial{trial:03}/cycle*/summary/cv_ranked.log | grep frame | awk '{{print $6}}'",
                shell=True,
                text=True,
                capture_output=True,
            )
            if res.returncode != 0:
                params.logger.error("error occurred at head command")
                sys.exit(1)
            y_list = [float(line) for line in res.stdout.strip().split("\n")]
            # y_list = []
            # print(f"trial{trial:03} is plotting")
            # for cycle in range(1, 210, 1):
            #     cycle_list = []
            #     for rep_path in Path(f"{params.feature_1d_directory}").glob(f"t{trial:03}c{cycle:03}r*.npy"):
            #         tmp_y_list = np.load(rep_path)
            #         cycle_list.append(max(tmp_y_list))
            #     else:
            #         if len(cycle_list) == 0:
            #             break
            #     y_list.append(max(cycle_list))
            plt.plot(y_list, label=f"trial{trial:03}", color=params.cmap(trial))
    plt.axhline(y=threshold, color="r", linestyle="--", label="unbound state")
    plt.legend()
    plt.xlabel("Cycle")
    plt.ylabel("Inter-COM [nm]")
    plt.grid()
    plt.tight_layout()
    # plt.xlim(0, 200)
    # plt.ylim(0, 10)
    plt.savefig(
        f"{params.output_directory}/images/cycle_feature.png",
        bbox_inches="tight",
        pad_inches=0.1,
        dpi=300,
    )
    if params.show_picture:
        plt.show()
    plt.clf()
    plt.close()


plot_feature(params=params, threshold=4)


# ## 1D: Clustering
# - Perform k-means clustering
def cluster_1d(
    params: Parameters,
    max_iter: int = 2000,
):
    for trial in params.n_trial_for_calc:
        # loads
        features = []
        a = np.load(f"{params.feature_1d_directory}/trial{trial}.npy")
        for each_rep in range(0, a.shape[0]):
            features_per_rep = a[each_rep]
            max_distance = np.max(features_per_rep)
            if max_distance > params.cutoff:
                continue
            else:
                features.append(features_per_rep)

        for n_clusters in params.n_clusters_for_try_1d:
            clustering_result_pkl = f"{params.output_directory}/cluster_objs/1d_trial{trial:03}_n_clusters{n_clusters}_cut{params.cutoff}.pkl"
            if Path(clustering_result_pkl).exists():
                params.logger.info(
                    f"Clustering result for trial{trial:03} n_clusters={n_clusters} cutoff={params.cutoff} already exists. Skipping."
                )
                continue

            # clustering
            params.logger.info(
                f"Starting clustering for trial{trial:03} n_clusters={n_clusters} cutoff={params.cutoff}"
            )
            estimator = deeptime.clustering.KMeans(
                n_clusters=n_clusters,
                max_iter=max_iter,
                metric="euclidean",
                tolerance=1e-05,
                init_strategy="kmeans++",
                fixed_seed=False,
                n_jobs=None,
                initial_centers=None,
                progress=None,
            )
            estimator.fit(np.concatenate(features))
            clustering_trjs = []
            for features_per_rep in features:
                clustering_trjs.append(estimator.transform(features_per_rep))
            clustering_model = estimator.fetch_model()

            # save result
            save_dict = {"trjs": clustering_trjs, "model": clustering_model}
            with open(clustering_result_pkl, "wb") as f:
                pickle.dump(save_dict, f)

            # create convergence figure
            inertias = clustering_model.inertias
            plt.figure(figsize=(5, 5))
            plt.plot(inertias)
            plt.xlabel("iteration")
            plt.ylabel("inertia")
            plt.title(f"trial{trial:03} n_clusters{n_clusters} clustering")
            plt.savefig(
                f"{params.output_directory}/images/clustering_converge_1d_trial{trial:03}_n_clusters{n_clusters}_cut{params.cutoff}.png",
                bbox_inches="tight",
                pad_inches=0.1,
                dpi=300,
            )
            if params.show_picture:
                plt.show()
            plt.clf()
            plt.close()

            params.logger.info(
                f"Finished clustering for trial{trial:03} n_clusters={n_clusters} cutoff={params.cutoff}"
            )


# run
cluster_1d(
    params=params,
)


# ## 1D: Plot Histgram of $d$
#
def plot_hist_1d_per_trial(
    params: Parameters,
    trial: int,
    n_clusters: int,
) -> None:
    # load clustering result
    params.logger.info(
        f"Starting histogram plotting for trial{trial:03} n_clusters{n_clusters}"
    )
    clustering_result_pkl = f"{params.output_directory}/cluster_objs/1d_trial{trial:03}_n_clusters{n_clusters}_cut{params.cutoff}.pkl"
    if not Path(clustering_result_pkl).exists():
        params.logger.info(
            f"Clustering obj for trial{trial:03} n_clusters={n_clusters} cutoff={params.cutoff} does not exist. Skipping. First, run clustering."
        )
        return
    with open(clustering_result_pkl, "rb") as f:
        save_dict = pickle.load(f)
        clustering_model = save_dict["model"]
        center_d = clustering_model.cluster_centers

    # load features
    features = []
    feature_file_trial = f"{params.feature_1d_directory}/trial{trial}.npy"
    a = np.load(feature_file_trial)
    for each_rep in range(0, a.shape[0]):
        features_per_rep = a[each_rep]
        max_distance = np.max(features_per_rep)
        if max_distance > params.cutoff:
            continue
        else:
            features.append(features_per_rep)

    # create figure
    plt.figure(figsize=(5, 3))
    plt.vlines(center_d, 0, 10**10, "tab:red", linewidth=0.5, label="cluster center")
    plt.vlines(  # comment out if you don't need cutoff line (default: cutoff=1000)
        params.cutoff, 0, 10**10, "magenta", linewidth=0.5, label="cutoff"
    )
    plt.hist(np.concatenate(features), bins=100, alpha=0.5)
    plt.yscale("log")
    plt.xlabel("$d$ [nm]")
    plt.ylabel("Frequency")
    plt.title(f"trial{trial:03} Inter-COM distance distribution")
    # put legend on right hand
    # plt.legend(loc='upper center')
    plt.ylim(1, 10**6)
    plt.savefig(
        f"{params.output_directory}/images/hist_1d_trial{trial:03}_n_clusters{n_clusters}_cut{params.cutoff}.png",
        bbox_inches="tight",
        pad_inches=0.1,
        dpi=300,
    )
    if params.show_picture:
        plt.show()
    plt.clf()
    plt.close()
    params.logger.info(
        f"Finished histogram plotting for trial{trial:03} n_clusters{n_clusters}"
    )


def plot_hist_1d(
    params: Parameters,
    n_clusters: int,
) -> None:
    for trial in params.n_trial_for_calc:
        plot_hist_1d_per_trial(
            params=params,
            trial=trial,
            n_clusters=n_clusters,
        )


# run
for n_clusters in params.n_clusters_for_try_1d:
    plot_hist_1d(params=params, n_clusters=n_clusters)


# ## 1D: Plot inertias
# - plot inertias along n_clusters
# - can be used to decide appropriate n_clusters
def plot_inertia_1d(params: Parameters):
    for trial in params.n_trial_for_calc:
        params.logger.info(f"Starting plotting inertia for trial{trial:03}")

        # load
        n_clusters_list = []
        inertias = []
        for n_clusters in params.n_clusters_for_try_1d:
            clustering_result = f"{params.output_directory}/cluster_objs/1d_trial{trial:03}_n_clusters{n_clusters}_cut{params.cutoff}.pkl"
            if not Path(clustering_result).exists():
                params.logger.info(
                    f"Clustering obj for trial{trial:03} n_clusters={n_clusters} cutoff={params.cutoff} does not exist. Skipping. First, run clustering."
                )
                continue
            with open(clustering_result, "rb") as f:
                save_dict = pickle.load(f)
                clustering_model = save_dict["model"]
            n_clusters_list.append(n_clusters)
            inertias.append(clustering_model.inertia)

        # check
        if len(n_clusters_list) == 0:
            params.logger.info(
                f"Clustering result not found for trial{trial:03}. No need to plot. Skipping."
            )
            return
        elif len(n_clusters_list) == 1:
            params.logger.info(
                f"Only one clustering result found for trial{trial:03}. Cannot plot. Skipping."
            )
            return

        # plot
        plt.figure(figsize=(3, 2))
        plt.plot(n_clusters_list, inertias, marker=".")
        plt.xlabel("number of clusters")
        plt.ylabel("inertia")
        plt.title(f"trial{trial:03}")
        plt.savefig(
            f"{params.output_directory}/images/inertia_1d_trial{trial:03}_cut{params.cutoff}.png",
            bbox_inches="tight",
            pad_inches=0.1,
            dpi=300,
        )
        if params.show_picture:
            plt.show()
        plt.clf()
        plt.close()
        params.logger.info(f"Finished plotting inertia for trial{trial:03}")


# run
plot_inertia_1d(params=params)


# ## 1D: Build MSM
# - build MSM for each trial and n_clusters as specified at the top.
def build_msm_1d(
    params: Parameters,
) -> None:
    for trial in params.n_trial_for_calc:
        for n_clusters in params.n_clusters_for_try_1d:
            # load or initialize msm result
            msm_result_pkl = f"{params.output_directory}/MSM_objs/1d_trial{trial:03}_n_clusters{n_clusters}_cut{params.cutoff}.pkl"
            if not Path(msm_result_pkl).exists():
                msm_result = dict()
            else:
                with open(msm_result_pkl, "rb") as f:
                    msm_result = pickle.load(f)

            # load or itinitialize count model result
            count_result_pkl = f"{params.output_directory}/Count_objs/1d_trial{trial:03}_n_clusters{n_clusters}_cut{params.cutoff}.pkl"
            if not Path(count_result_pkl).exists():
                count_result = dict()
            else:
                with open(count_result_pkl, "rb") as f:
                    count_result = pickle.load(f)

            # check if all clustering results exist
            clustering_result_pkl = f"{params.output_directory}/cluster_objs/1d_trial{trial:03}_n_clusters{n_clusters}_cut{params.cutoff}.pkl"
            if not Path(clustering_result_pkl).exists():
                params.logger.info(
                    f"Clustering obj for trial{trial:03} n_clusters={n_clusters} cutoff={params.cutoff} does not exist. Skipping. First, run clustering."
                )
                continue
            with open(clustering_result_pkl, "rb") as f:
                save_dict = pickle.load(f)
                trjs = save_dict["trjs"]
                clustering_model = save_dict["model"]

            # build msm for all lags in params.lags_for_try_1d
            for lag in params.lags_for_try_1d:
                # skip if already exits
                if lag in msm_result.keys():
                    params.logger.info(
                        f"MSM result for trial{trial:03} n_clusters={n_clusters} lag={lag} already exists. Skipping."
                    )
                    continue

                params.logger.info(
                    f"Starting MSM for trial{trial:03} n_clusters={n_clusters} lag={lag}"
                )
                # build count matrix
                count_estimator = deeptime.markov.TransitionCountEstimator(
                    lagtime=lag, count_mode="sliding", n_states=None, sparse=False
                )
                count_model = count_estimator.fit(trjs).fetch_model()

                # build msm
                msm_estimator = deeptime.markov.msm.MaximumLikelihoodMSM(
                    reversible=True,
                    stationary_distribution_constraint=None,
                    sparse=False,
                    allow_disconnected=False,
                    maxiter=1000000,
                    maxerr=1e-08,
                    connectivity_threshold=0,
                    transition_matrix_tolerance=1e-06,
                    lagtime=None,
                    use_lcc=False,
                )
                try:
                    msm_model = msm_estimator.fit(count_model).fetch_model()
                except Exception as e:
                    params.logger.error(
                        f"Error occurred in building MSM at trial{trial:03} n_clusters={n_clusters} lag={lag}"
                    )
                    params.logger.error(f"error message: {e}")
                    continue

                # save msm result
                msm_result[lag] = msm_model
                count_result[lag] = count_model

                with open(msm_result_pkl, "wb") as f:
                    pickle.dump(msm_result, f)
                with open(count_result_pkl, "wb") as f:
                    pickle.dump(count_result, f)

                # save cluster center coordinates and corresponding stationary destribution
                n_clusters_in_clustering = len(clustering_model.cluster_centers)
                n_clusters_in_msm = msm_model.n_states
                cluster_centers = clustering_model.cluster_centers

                stationary_distribution = np.zeros(
                    n_clusters_in_clustering,
                )
                largest_connected_set = (
                    deeptime.markov.tools.estimation.largest_connected_set(
                        count_model.count_matrix, directed=True
                    )
                )
                stationary_distribution[largest_connected_set] = (
                    msm_model.stationary_distribution
                )

                if n_clusters_in_clustering != n_clusters_in_msm:
                    params.logger.info(
                        f"The number of clusters for clustering and msm are different for trial{trial:03} n_clusters={n_clusters} lag={lag}"
                    )
                    params.logger.info(
                        f"Number of clusters for clustering: {n_clusters_in_clustering}"
                    )
                    params.logger.info(
                        f"Number of clusters in the largest connected set: {n_clusters_in_msm}"
                    )

                center_pi_csv = f"{params.output_directory}/result_csvs/1d_trial{trial:03}_n_clusters{n_clusters}_lag{lag}_cut{params.cutoff}.csv"
                with open(center_pi_csv, "w") as f:
                    f.write("cluster_center_d,cluster_center_pi\n")
                    num_clusters = len(cluster_centers)
                    for cluster_ind in range(num_clusters):
                        center_d = cluster_centers[cluster_ind][0]
                        pi = stationary_distribution[cluster_ind]
                        f.write(f"{center_d},{pi}\n")
                params.logger.info(
                    f"Finished MSM for trial{trial:03} n_clusters={n_clusters} lag={lag}"
                )
            params.logger.info(
                f"Finished MSM for trial{trial:03} n_clusters={n_clusters}"
            )


# run
build_msm_1d(params=params)


# ## 1D: Plot ITS using distance
# - calculate "Implied Time Scale" for each trial
def plot_its_1d(params: Parameters) -> None:
    for trial in params.n_trial_for_calc:
        for n_clusters in params.n_clusters_for_try_1d:
            params.logger.info(
                f"Starting plotting ITS for trial{trial:03} n_clusters={n_clusters}"
            )
            # load or initialize msm result
            msm_result_pkl = f"{params.output_directory}/MSM_objs/1d_trial{trial:03}_n_clusters{n_clusters}_cut{params.cutoff}.pkl"
            if not Path(msm_result_pkl).exists():
                msm_result = dict()
            else:
                with open(msm_result_pkl, "rb") as f:
                    msm_result = pickle.load(f)

            # accumulate msm models and calculate implied timescales
            msm_models = []
            for lag in sorted(msm_result.keys()):
                msm_models.append(msm_result[lag])
            its_data = deeptime.util.validation.implied_timescales(
                models=msm_models,
                n_its=10,  # decrease this number for n_clusters < 11
            )

            # plot its
            plt.figure(figsize=(3, 2))
            deeptime.plots.plot_implied_timescales(
                its_data,
                n_its=10,  # decrease this number for n_clusters < 11
                process=None,
                show_mle=True,
                show_sample_mean=True,  # ignored
                show_sample_confidence=True,  # ignored
                show_cutoff=True,  # ignored
                sample_confidence=0.95,  # ignored
                colors=None,
                # ax=None,
            )
            plt.title(f"trial{trial:03} n_clusters{n_clusters}")
            plt.xlabel("lag time (steps)")  # 1 step = 1 interval in the trajectory
            plt.ylabel(
                "implied timescale (steps)"
            )  # 1 step = 1 interval in the trajectory
            plt.yscale("log")
            plt.savefig(
                f"{params.output_directory}/images/its_1d_trial{trial:03}_n_clusters{n_clusters}_cut{params.cutoff}.png",
                bbox_inches="tight",
                pad_inches=0.1,
                dpi=300,
            )
            if params.show_picture:
                plt.show()
            plt.clf()
            plt.close()
            params.logger.info(
                f"Finished plotting ITS for trial{trial:03} n_clusters={n_clusters}"
            )


# run
plot_its_1d(params=params)


def plot_fel_along_d_1d(
    params: Parameters,
    n_clusters: int,
    lag: int,
    interpolate_type: str = "linear",
) -> None:
    # check interpolation type
    if interpolate_type not in ["linear", "cubic", "none"]:
        params.logger.error(f"Interpolation type {interpolate_type} is not supported.")
        return

    # log
    params.logger.info(
        f"Starting plotting FEL along d for n_clusters={n_clusters} lag={lag}"
    )

    if interpolate_type == "none":
        # doesn't interpolate for each trial
        plot_data_x_all = []
        plot_data_y_all = []
        plt.figure(figsize=(5, 5))
        for trial in params.n_trial_for_calc:
            print(f"{trial=}")
            print(f"{n_clusters=}")
            # load
            msm_result_csv = f"{params.output_directory}/result_csvs/1d_trial{trial:03}_n_clusters{n_clusters}_lag{lag}_cut{params.cutoff}.csv"
            if not Path(msm_result_csv).exists():
                params.logger.info(
                    "MSM result csv file not found. Skipping. First, run MSM."
                )
                continue
            df = pd.read_csv(msm_result_csv)
            cluster_center_d = df["cluster_center_d"].values
            stat_dists = df["cluster_center_pi"].values

            # select clusters where stationary distribution is not zero
            mask = stat_dists > 0
            cluster_center_d = cluster_center_d[mask]
            stat_dists = stat_dists[mask]

            # sort by d
            x_args = np.argsort(cluster_center_d)
            cluster_center_d = cluster_center_d[x_args]
            stat_dists = stat_dists[x_args]

            # plot
            y_values = params._RT * np.log(stat_dists)
            y_values -= np.min(y_values)
            plt.plot(
                cluster_center_d,
                y_values,
                label=f"trial{trial:03}",
                color=params.cmap(trial),
            )
            plot_data_x_all.extend(cluster_center_d.tolist())
            plot_data_y_all.extend(y_values.tolist())

        # define bins
        plot_data_x_all = np.array(plot_data_x_all)
        plot_data_y_all = np.array(plot_data_y_all)
        print(f"{plot_data_x_all=}")
        print(f"{plot_data_y_all=}")
        if len(plot_data_x_all) == 0:
            params.logger.info("No data to plot. Skipping.")
            return
        bins = np.linspace(
            start=min(plot_data_x_all), stop=max(plot_data_x_all), num=params.nbins
        )
        bins -= (bins[1] - bins[0]) / 2
        bin_centers = bins[:-1] + (bins[1] - bins[0]) / 2
        bin_stat_dist_means = np.zeros(len(bins) - 1)
        bin_stat_dist_stds = np.zeros(len(bins) - 1)
        for i in range(len(bins) - 1):
            bin_mask = (plot_data_x_all >= bins[i]) & (plot_data_x_all < bins[i + 1])
            bin_stat_dist_means[i] = np.mean(plot_data_y_all[bin_mask])
            bin_stat_dist_stds[i] = np.std(plot_data_y_all[bin_mask])
        bin_stat_dist_stes = bin_stat_dist_stds / np.sqrt(len(params.n_trial_for_calc))

        # bin_stat_dist_means -= np.min(bin_stat_dist_means)

        plt.plot(
            bin_centers,
            bin_stat_dist_means,
            color="black",
            linewidth=2,
        )
        plt.errorbar(
            bin_centers,
            bin_stat_dist_means,
            yerr=bin_stat_dist_stes,
            fmt="o",
            color="black",
            label=r"Mean $\pm$SE",
            capsize=3,
        )

        plt.xlim(0, params.cutoff)
        plt.ylim(0, 20)
        plt.title(f"n_clusters={n_clusters} lag={lag}")
        plt.xlabel(r"$d$ [nm]")
        plt.ylabel(r"$\Delta W$ [kcal/mol]")
        # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig(
            f"{params.output_directory}/images/fel_1d_n_clusters{n_clusters}_lag{lag}_interp_{interpolate_type}.png",
            bbox_inches="tight",
            pad_inches=0.1,
            dpi=300,
        )
        if params.show_picture:
            plt.show()
        plt.clf()
        plt.close()

    elif interpolate_type in ["linear", "cubic"]:
        # interpolate for each trial
        plot_data_x_all = []
        plot_data_y_all = []
        plt.figure(figsize=(5, 5))

        # define bins
        center_list = []
        for trial in params.n_trial_for_calc:
            # load
            msm_result_csv = f"{params.output_directory}/result_csvs/1d_trial{trial:03}_n_clusters{n_clusters}_lag{lag}_cut{params.cutoff}.csv"
            if not Path(msm_result_csv).exists():
                params.logger.info(
                    "MSM result csv file not found. Skipping. First, run MSM."
                )
                continue
            df = pd.read_csv(msm_result_csv)
            cluster_center_d = df["cluster_center_d"].values.tolist()
            center_list.extend(cluster_center_d)
        if len(center_list) == 0:
            params.logger.info("No data to plot. Skipping.")
            return
        bins = np.linspace(
            start=min(center_list), stop=max(center_list), num=params.nbins
        )
        bins -= (bins[1] - bins[0]) / 2
        bin_centers = bins[:-1] + (bins[1] - bins[0]) / 2
        n_intervals = len(bins) - 1
        n_trials = len(params.n_trial_for_calc)
        bin_stat_dist_trials = np.zeros((n_intervals, n_trials))
        y_values_trials = np.zeros((n_intervals, n_trials))

        for trial_index, trial in enumerate(params.n_trial_for_calc):
            # load
            msm_result_csv = f"{params.output_directory}/result_csvs/1d_trial{trial:03}_n_clusters{n_clusters}_lag{lag}_cut{params.cutoff}.csv"
            if not Path(msm_result_csv).exists():
                params.logger.info(
                    "MSM result csv file not found. Skipping. First, run MSM."
                )
                continue
            df = pd.read_csv(msm_result_csv)
            cluster_center_d = df["cluster_center_d"].values
            stat_dists = df["cluster_center_pi"].values

            # summation for all intervals
            for bin_i in range(n_intervals):
                mask = (cluster_center_d >= bins[bin_i]) & (
                    cluster_center_d < bins[bin_i + 1]
                )
                bin_stat_dist_trials[bin_i, trial_index] = np.sum(stat_dists[mask])

            # interpolate for the missing values
            zero_mask = bin_stat_dist_trials[:, trial_index] == 0
            nonzero_mask = bin_stat_dist_trials[:, trial_index] > 0

            x_values_nonzero = bin_centers[nonzero_mask]
            y_values_nonzero = params._RT * np.log(
                bin_stat_dist_trials[nonzero_mask, trial_index]
            )

            x_values_zero = bin_centers[zero_mask]

            if interpolate_type == "cubic":
                cs = CubicSpline(
                    x_values_nonzero,
                    y_values_nonzero,
                    bc_type="natural",
                    extrapolate=True,
                )
                y_values_zero = cs(x_values_zero)
            elif interpolate_type == "linear":
                y_values_zero = np.interp(
                    x_values_zero,
                    x_values_nonzero,
                    y_values_nonzero,
                )

            n_points_to_interpolate = len(y_values_zero)
            params.logger.info(
                f"trial{trial:03} interpolated {n_points_to_interpolate} empty bin(s)"
                f" out of {len(bin_centers)} bins by {interpolate_type} interpolation"
            )

            y_values_trials[zero_mask, trial_index] = y_values_zero
            y_values_trials[nonzero_mask, trial_index] = y_values_nonzero

            # plot
            y_values_trials[:, trial_index] -= np.min(y_values_trials[:, trial_index])
            # y_values = y_values_trials[:, trial_index].copy()
            # y_values -= np.min(y_values)
            plt.plot(
                bin_centers,
                y_values_trials[:, trial_index],
                label=f"trial{trial:03}",
                color=params.cmap(trial),
            )

        # skip MSM-failing trials to calculate mean and std
        msm_ok_trial_indices = np.where(~np.all(y_values_trials == 0, axis=0))[0]
        y_values_trials = y_values_trials[:, msm_ok_trial_indices]
        n_trials_ok = len(msm_ok_trial_indices)

        # statistics
        bin_stat_dist_means = np.mean(y_values_trials, axis=1)
        bin_stat_dist_stds = np.std(y_values_trials, axis=1)
        bin_stat_dist_stes = bin_stat_dist_stds / np.sqrt(n_trials_ok)

        plt.plot(
            bin_centers,
            bin_stat_dist_means,
            color="black",
            linewidth=2,
        )
        plt.errorbar(
            bin_centers,
            bin_stat_dist_means,
            yerr=bin_stat_dist_stes,
            fmt="o",
            color="black",
            label=r"Mean $\pm$SE",
            capsize=3,
        )

        plt.xlim(0, params.cutoff)
        plt.ylim(0, 20)
        plt.title(f"n_clusters={n_clusters} lag={lag}")
        plt.xlabel(r"$d$ [nm]")
        plt.ylabel(r"$\Delta W$ [kcal/mol]")
        # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig(
            f"{params.output_directory}/images/fel_1d_n_clusters{n_clusters}_lag{lag}_interp_{interpolate_type}.png",
            bbox_inches="tight",
            pad_inches=0.1,
            dpi=300,
        )
        if params.show_picture:
            plt.show()
        plt.clf()
        plt.close()

    else:
        params.logger.error(f"Interpolation type {interpolate_type} is not supported.")
        return

    params.logger.info(
        f"Finished plotting FEL along d for n_clusters={n_clusters} lag={lag}"
    )


# run
for n_clusters in n_clusters_1d_main:
    for lag in lags_1d_main:
        plot_fel_along_d_1d(
            params=params, n_clusters=n_clusters, lag=lag, interpolate_type="none"
        )
        plot_fel_along_d_1d(
            params=params,
            n_clusters=n_clusters,
            lag=lag,
        )
        plot_fel_along_d_1d(
            params=params, n_clusters=n_clusters, lag=lag, interpolate_type="cubic"
        )


def calc_vc_per_trial(
    params: Parameters,
    trial: int,
    lower_unbound: float,
    upper_unbound: float,
):
    # check the order of bound, unbound
    if lower_unbound > upper_unbound:
        params.logger.error("lower_unbound must be less than upper_unbound. Exiting.")
        sys.exit(1)
    if upper_unbound > params.cutoff:
        params.logger.error("upper_unbound must be less than cut_off. Exiting.")
        sys.exit(1)

    # check if csv file already exists
    vc_csv = f"{params.output_directory}/result_csvs/1d_VC.csv"
    if not Path(vc_csv).exists():
        with open(vc_csv, "w") as f:
            f.write("trial,lower_unbound,upper_unbound,volume,correction\n")
    else:
        df = pd.read_csv(vc_csv)
        if (
            len(
                df[
                    (df["trial"] == trial)
                    & (df["lower_unbound"] == lower_unbound)
                    & (df["upper_unbound"] == upper_unbound)
                ]
            )
            > 0
        ):
            # debug level to avoid too noisy output
            params.logger.debug(
                f"Volume correction for trial{trial:03} lower_unbound{lower_unbound} upper_unbound{upper_unbound} already calculated. Skipping"
            )
            return

    params.logger.info(f"Starting volume correction calculation for trial{trial:03}")
    # load
    inter_COM_vec_all = []
    feature_file_trial = f"{params.feature_3d_directory}/trial{trial}.npy"
    a = np.load(feature_file_trial)
    for each_rep in range(0, a.shape[0]):
        inter_COM_vec = a[each_rep]
        max_distance = np.max(np.linalg.norm(inter_COM_vec, axis=1))
        if max_distance > params.cutoff:
            continue
        inter_COM_vec_all.append(inter_COM_vec)
    inter_COM_vec_all = np.concatenate(inter_COM_vec_all)

    # select snapshots of unbound states
    inter_COM_d_all = np.linalg.norm(inter_COM_vec_all, axis=1)
    ix = np.where(
        (inter_COM_d_all >= lower_unbound) & (inter_COM_d_all < upper_unbound)
    )
    inter_COM_vec = inter_COM_vec_all[ix]

    # calculate volume correction
    if len(inter_COM_vec) > 0:
        hull = ConvexHull(inter_COM_vec)
        volume = hull.volume * 1000  # [angstrome^3]
        correction = params._RT * np.log(volume / 1661)  # [kcal/mol]
    else:
        volume = 0
        correction = 0

    # save to csv file
    with open(vc_csv, "a") as f:
        f.write(f"{trial},{lower_unbound},{upper_unbound},{volume},{correction}\n")

    params.logger.info(f"Finished volume correction calculation for trial{trial:03}")


def calc_binding_energy_1d(
    params: Parameters,
    n_clusters: int,
    lag: int,
    lower_bound: float,
    upper_bound: float,
    lower_unbound: float,
    upper_unbound: float,
) -> None:
    # check the order of bound, unbound
    if lower_bound >= upper_bound:
        params.logger.error("lower_bound should be less than upper_bound. Exiting.")
        sys.exit(1)
    if lower_unbound >= upper_unbound:
        params.logger.error("lower_unbound should be less than upper_unbound. Exiting.")
        sys.exit(1)
    if upper_bound > lower_unbound:
        params.logger.error(
            "upper_bound should be less than or equal to lower_unbound. Exiting."
        )
        sys.exit(1)
    if upper_unbound > params.cutoff:
        params.logger.error(
            "upper_unbound should be less than or equal to cutoff. Exiting."
        )
        sys.exit(1)

    # check if csv file already exists
    be_csv = f"{params.output_directory}/result_csvs/1d_binding_energy_summary.csv"
    if not Path(be_csv).exists():
        with open(be_csv, "w") as f:
            f.write(
                "trial,n_clusters,lag,lower_bound,upper_bound,lower_unbound,upper_unbound,dG_PMF,VC\n"
            )

    # calc volume correction
    if params.do_volume_correction:
        for trial in params.n_trial_for_calc:
            calc_vc_per_trial(
                params=params,
                trial=trial,
                lower_unbound=lower_unbound,
                upper_unbound=upper_unbound,
            )

    # calc binding energy
    for trial in params.n_trial_for_calc:
        # check if already calculated
        df = pd.read_csv(be_csv)
        if (
            len(
                df[
                    (df["trial"] == trial)
                    & (df["n_clusters"] == n_clusters)
                    & (df["lag"] == lag)
                    & (df["lower_bound"] == lower_bound)
                    & (df["upper_bound"] == upper_bound)
                    & (df["lower_unbound"] == lower_unbound)
                    & (df["upper_unbound"] == upper_unbound)
                ]
            )
            > 0
        ):
            params.logger.debug(
                f"Binding energy for trial{trial:03} "
                f"n_clusters{n_clusters} "
                f"lag{lag} "
                f"lower_bound{lower_bound} "
                f"upper_bound{upper_bound} "
                f"lower_unbound{lower_unbound} "
                f"upper_unbound{upper_unbound} "
                "already exists. Skipped"
            )
            continue

        # load msm
        msm_result_csv = f"{params.output_directory}/result_csvs/1d_trial{trial:03}_n_clusters{n_clusters}_lag{lag}_cut{params.cutoff}.csv"
        if not Path(msm_result_csv).exists():
            params.logger.info(
                "MSM result csv file not found. Skipping. First, run MSM."
            )
            continue
        stat_dist_df = pd.read_csv(msm_result_csv)
        cluster_center_d = stat_dist_df["cluster_center_d"].values
        stat_dists = stat_dist_df["cluster_center_pi"].values

        # select clusters where stationary distribution is not zero
        mask = stat_dists > 0
        cluster_center_d = cluster_center_d[mask]
        stat_dists = stat_dists[mask]

        # select clusters in bound/unbound states
        ix_b = np.where(
            (cluster_center_d >= lower_bound) & (cluster_center_d < upper_bound)
        )
        ix_u = np.where(
            (cluster_center_d >= lower_unbound) & (cluster_center_d < upper_unbound)
        )
        pi_b = np.sum(stat_dists[ix_b])
        pi_u = np.sum(stat_dists[ix_u])

        # check if no cluster found in either bound/unbound state
        if len(stat_dists[ix_b]) == 0 or len(stat_dists[ix_u]) == 0:
            params.logger.info(
                "No cluster found in either bound/unbound state for "
                f"trial{trial:03} "
                f"n_clusters={n_clusters} "
                f"lag={lag} "
                f"lower_bound={lower_bound} "
                f"upper_bound={upper_bound} "
                f"lower_unbound={lower_unbound} "
                f"upper_unbound={upper_unbound} . Skipping."
                "Please adjust the bound/unbound state definition"
            )
            continue
        dG = params._RT * np.log(pi_b / pi_u)  # [kcal/mol]

        # get vc information
        if params.do_volume_correction:
            vc_csv = f"{params.output_directory}/result_csvs/1d_VC.csv"
            vc_df = pd.read_csv(vc_csv)
            vc_df = vc_df[
                (vc_df["lower_unbound"] == lower_unbound)
                & (vc_df["upper_unbound"] == upper_unbound)
            ]
            vc = vc_df[vc_df["trial"] == trial]["correction"].values[0]
        else:
            vc = 0

        # save to csv file
        with open(be_csv, "a") as f:
            f.write(
                f"{trial},{n_clusters},{lag},{lower_bound},{upper_bound},{lower_unbound},{upper_unbound},{dG},{vc}\n"
            )

        # debug level to avoid too noisy output
        params.logger.debug(
            f"Finished calculating binding energy for "
            f"trial{trial:03} "
            f"n_clusters={n_clusters} "
            f"lag={lag} "
            f"lower_bound={lower_bound} "
            f"upper_bound={upper_bound} "
            f"lower_unbound={lower_unbound} "
            f"upper_unbound={upper_unbound} "
        )
    params.logger.info(
        f"Finished calculating binding energy for all trials for "
        f"n_clusters={n_clusters} "
        f"lag={lag} "
        f"lower_bound={lower_bound} "
        f"upper_bound={upper_bound} "
        f"lower_unbound={lower_unbound} "
        f"upper_unbound={upper_unbound} "
    )


# calc binding energy
for n_clusters in n_clusters_1d_main:
    for lag in lags_1d_main:
        calc_binding_energy_1d(
            params=params,
            n_clusters=n_clusters,
            lag=lag,
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            lower_unbound=lower_unbound,
            upper_unbound=params.cutoff,
        )


# 3D: Clustering
def cluster_3d(
    params: Parameters,
    max_iter: int = 2000,
):
    for n_clusters in params.n_clusters_for_try_3d:
        # check
        clustering_result_pkl = f"{params.output_directory}/cluster_objs/3d_n_clusters{n_clusters}_cut{params.cutoff}.pkl"
        if Path(clustering_result_pkl).exists():
            params.logger.info(
                f"Clustering result obj for n_clusters={n_clusters} cutoff={params.cutoff} already exists. Skipping."
            )
            continue

        # cluster
        params.logger.info(
            f"Starting clustering for n_clusters={n_clusters} cutoff={params.cutoff}"
        )
        # loads
        features = []
        for trial in params.n_trial_for_calc:
            a = np.load(f"{params.feature_3d_directory}/trial{trial}.npy")
            for each_rep in range(0, a.shape[0]):
                feature_per_rep = a[each_rep]
                max_distance = np.max(np.linalg.norm(feature_per_rep, axis=1))
                if max_distance > params.cutoff:
                    continue
                else:
                    features.append(feature_per_rep)

        # clustering
        estimator = deeptime.clustering.KMeans(
            n_clusters=n_clusters,
            max_iter=max_iter,
            metric="euclidean",
            tolerance=1e-05,
            init_strategy="kmeans++",
            fixed_seed=False,
            n_jobs=None,
            initial_centers=None,
            progress=None,
        )
        estimator.fit(np.concatenate(features))
        clustering_trjs = []
        for features_per_rep in features:
            clustering_trjs.append(estimator.transform(features_per_rep))
        clustering_model = estimator.fetch_model()

        # save result
        save_dict = {"trjs": clustering_trjs, "model": clustering_model}
        with open(clustering_result_pkl, "wb") as f:
            pickle.dump(save_dict, f)

        # create convergence figure
        inertias = clustering_model.inertias
        plt.figure(figsize=(5, 5))
        plt.plot(inertias)
        plt.xlabel("iteration")
        plt.ylabel("inertia")
        plt.title(f"n_clusters{n_clusters} clustering")
        plt.savefig(
            f"{params.output_directory}/images/clustering_converge_3d_n_clusters{n_clusters}_cut{params.cutoff}.png",
            bbox_inches="tight",
            pad_inches=0.1,
            dpi=300,
        )
        if params.show_picture:
            plt.show()
        plt.clf()
        plt.close()

        params.logger.info(
            f"Finished clustering for n_clusters={n_clusters} cutoff={params.cutoff}"
        )


# plot
cluster_3d(params=params)


# 3D: Plot Histgram of $d$
def plot_hist_3d(
    params: Parameters,
    n_clusters: int,
) -> None:
    # check if clustering result exists
    clustering_result_pkl = f"{params.output_directory}/cluster_objs/3d_n_clusters{n_clusters}_cut{params.cutoff}.pkl"
    if not Path(clustering_result_pkl).exists():
        params.logger.info(
            f"Clustering obj for n_clusters={n_clusters} cutoff={params.cutoff} does not exist. Skipping. First, run clustering."
        )
        return
    with open(clustering_result_pkl, "rb") as f:
        save_dict = pickle.load(f)
        clustering_model = save_dict["model"]
        center_d = np.linalg.norm(clustering_model.cluster_centers, axis=1)

    # load
    features = []
    for trial in params.n_trial_for_calc:
        a = np.load(f"{params.feature_3d_directory}/trial{trial}.npy")
        for each_rep in range(0, a.shape[0]):
            features_per_rep = a[each_rep]
            max_distance = np.max(np.linalg.norm(features_per_rep, axis=1))
            if max_distance > params.cutoff:
                continue
            else:
                features.append(np.linalg.norm(features_per_rep, axis=1))

    # create figure
    plt.figure(figsize=(5, 3))
    plt.vlines(center_d, 0, 10**10, "tab:red", linewidth=0.5, label="cluster center")
    plt.vlines(params.cutoff, 0, 10**10, "magenta", linewidth=0.5, label="cutoff")
    plt.hist(np.concatenate(features), bins=100, alpha=0.5)
    plt.yscale("log")
    plt.xlabel("$d$ [nm]")
    plt.ylabel("Frequency")
    plt.title("Inter-COM distance distribution")
    # put legend on right hand
    plt.legend(loc="upper center")
    plt.ylim(1, 10**6)
    plt.savefig(
        f"{params.output_directory}/images/hist_3d_n_clusters{n_clusters}_cut{params.cutoff}.png",
        bbox_inches="tight",
        pad_inches=0.1,
        dpi=300,
    )
    if params.show_picture:
        plt.show()
    plt.clf()
    plt.close()

    params.logger.info(f"Finished histogram plotting for n_clusters={n_clusters}")


# run
for n_clusters in params.n_clusters_for_try_3d:
    plot_hist_3d(params=params, n_clusters=n_clusters)


# 3D: Plot inertias
# - plot inertias along n_clusters
# - can be used to decide appropriate n_clusters
def plot_inertia_3d(params: Parameters):
    n_clusters_list = []
    inertias = []
    for n_clusters in params.n_clusters_for_try_3d:
        clustering_result_pkl = f"{params.output_directory}/cluster_objs/3d_n_clusters{n_clusters}_cut{params.cutoff}.pkl"
        if not Path(clustering_result_pkl).exists():
            params.logger.info(
                f"Clustering result obj for n_clusters={n_clusters} cutoff={params.cutoff} does not exist. Skipping."
            )
            continue
        with open(clustering_result_pkl, "rb") as f:
            save_dict = pickle.load(f)
            clustering_model = save_dict["model"]
        n_clusters_list.append(n_clusters)
        inertias.append(clustering_model.inertia)

    # check
    if len(n_clusters_list) == 0:
        params.logger.info(
            "No clustering result found for 3d. No need to plot. Skipping."
        )
        return
    elif len(n_clusters_list) == 1:
        params.logger.info("Only one clustering result found. Cannot plot. Skipping.")
        return

    # plot
    params.logger.info("Starting plotting inertia for 3d")
    plt.figure(figsize=(3, 2))
    plt.plot(n_clusters_list, inertias, marker=".")
    plt.xlabel("number of clusters")
    plt.ylabel("inertia")
    plt.title("3d")
    plt.savefig(
        f"{params.output_directory}/images/inertia_3d_cut{params.cutoff}.png",
        bbox_inches="tight",
        pad_inches=0.1,
        dpi=300,
    )
    if params.show_picture:
        plt.show()
    plt.clf()
    plt.close()

    params.logger.info("Finished plotting inertia for 3d")


# run
plot_inertia_3d(params=params)


# 3D: Build MSM
def build_msm_3d(
    params: Parameters,
) -> None:
    for n_clusters in params.n_clusters_for_try_3d:
        # load or initialize msm result
        msm_result_pkl = f"{params.output_directory}/MSM_objs/3d_n_clusters{n_clusters}_cut{params.cutoff}.pkl"
        if not Path(msm_result_pkl).exists():
            msm_result = dict()
        else:
            with open(msm_result_pkl, "rb") as f:
                msm_result = pickle.load(f)

        # load or itinitialize count model result
        count_result_pkl = f"{params.output_directory}/Count_objs/3d_n_clusters{n_clusters}_cut{params.cutoff}.pkl"
        if not Path(count_result_pkl).exists():
            count_result = dict()
        else:
            with open(count_result_pkl, "rb") as f:
                count_result = pickle.load(f)

        # check if all clustering results exist
        clustering_result_pkl = f"{params.output_directory}/cluster_objs/3d_n_clusters{n_clusters}_cut{params.cutoff}.pkl"
        if not Path(clustering_result_pkl).exists():
            params.logger.info(
                f"Clustering obj for n_clusters={n_clusters} cutoff={params.cutoff} does not exist. Skipping. First, run clustering."
            )
            continue
        with open(clustering_result_pkl, "rb") as f:
            save_dict = pickle.load(f)
            trjs = save_dict["trjs"]
            clustering_model = save_dict["model"]

        # build msm for all lags in params.lags_for_try_3d
        for lag in params.lags_for_try_3d:
            # skip if already exits
            if lag in msm_result.keys():
                params.logger.info(
                    f"MSM result for n_clusters={n_clusters} lag={lag} already exists. Skipping."
                )
                continue

            params.logger.info(f"Starting MSM for n_clusters={n_clusters} lag={lag}")
            # build count matrix
            count_estimator = deeptime.markov.TransitionCountEstimator(
                lagtime=lag, count_mode="effective", n_states=None, sparse=False
            )
            count_model = count_estimator.fit(
                trjs,
                # "n_jobs = None" means using all cores.
                # If you want to save the memory at the cost of speed,
                # specify n_jobs = integer meaning maximum number ofcores to use.
                n_jobs=None,
            ).fetch_model()

            # build msm
            msm_estimator = deeptime.markov.msm.BayesianMSM(
                n_samples=100,
                n_steps=None,  # 1 step = 1 interval in the trajectory
                reversible=True,
                stationary_distribution_constraint=None,
                sparse=False,
                maxiter=1000000,
                maxerr=1e-08,
                lagtime=None,
            )
            try:
                msm_model = msm_estimator.fit(count_model).fetch_model()
            except Exception as e:
                params.logger.error(
                    f"Error occurred in building MSM at n_clusters={n_clusters} lag={lag}"
                )
                params.logger.error(f"error message: {e}")
                continue

            # save msm result
            msm_result[lag] = msm_model
            count_result[lag] = count_model
            with open(msm_result_pkl, "wb") as f:
                pickle.dump(msm_result, f)
            with open(count_result_pkl, "wb") as f:
                pickle.dump(count_result, f)

            # save cluster center coordinates and corresponding stationary destribution
            n_clusters_in_clustering = len(clustering_model.cluster_centers)
            n_clusters_in_msm = msm_model.prior.n_states
            cluster_centers = clustering_model.cluster_centers

            stationary_distribution = np.zeros(
                (len(msm_model.samples), n_clusters_in_clustering)
            )
            largest_connected_set = (
                deeptime.markov.tools.estimation.largest_connected_set(
                    count_model.count_matrix, directed=True
                )
            )
            for i, sample in enumerate(msm_model.samples):
                stationary_distribution[i, largest_connected_set] = (
                    sample.stationary_distribution
                )
            if n_clusters_in_clustering != n_clusters_in_msm:
                params.logger.info(
                    f"The number of clusters for clustering and msm are different for n_clusters={n_clusters} lag={lag}"
                )
                params.logger.info(
                    f"Number of clusters for clustering: {n_clusters_in_clustering}"
                )
                params.logger.info(
                    f"Number of clusters in the largest connected set: {n_clusters_in_msm}"
                )
            center_pi_csv = f"{params.output_directory}/result_csvs/3d_n_clusters{n_clusters}_lag{lag}_cut{params.cutoff}.csv"
            with open(center_pi_csv, "w") as f:
                header = "cluster_center_x,cluster_center_y,cluster_center_z,"
                for i in range(len(msm_model.samples)):
                    header += f"cluster_center_pi_{i + 1},"
                header = header[:-1]
                header += "\n"
                f.write(header)
                num_clusters = len(cluster_centers)
                for cluster_ind in range(num_clusters):
                    center_x = cluster_centers[cluster_ind][0]
                    center_y = cluster_centers[cluster_ind][1]
                    center_z = cluster_centers[cluster_ind][2]
                    pis_in_samples = stationary_distribution[:, cluster_ind]
                    pis_str = ",".join([str(pi) for pi in pis_in_samples])[:-1]
                    f.write(f"{center_x},{center_y},{center_z},{pis_str}\n")
            params.logger.info(f"Finished MSM for n_clusters={n_clusters} lag={lag}")


# build msm
build_msm_3d(params)


# 3D: Plot ITS
def plot_its_3d(
    params: Parameters,
) -> None:
    for n_clusters in params.n_clusters_for_try_3d:
        # load or initialize msm result
        msm_result_pkl = f"{params.output_directory}/MSM_objs/3d_n_clusters{n_clusters}_cut{params.cutoff}.pkl"
        if not Path(msm_result_pkl).exists():
            params.logger.info(
                f"MSM result obj for n_clusters={n_clusters} cutoff={params.cutoff} does not exist. Skipping."
                f"First, run MSM."
            )
            continue
        else:
            with open(msm_result_pkl, "rb") as f:
                msm_result = pickle.load(f)

        # accumulate msm models and calculate implied timescales
        params.logger.info(f"Starting plotting ITS for n_clusters={n_clusters}")
        msm_models = []
        for lag in sorted(msm_result.keys()):
            msm_models.append(msm_result[lag])
        its_data = deeptime.util.validation.implied_timescales(
            models=msm_models,
            n_its=10,  # decrease this number for n_clusters < 11
        )

        # plot its
        plt.figure(figsize=(3, 2))
        deeptime.plots.plot_implied_timescales(
            its_data,
            n_its=10,  # decrease this number for n_clusters < 11
            process=None,
            show_mle=True,
            show_sample_mean=True,  # ignored
            show_sample_confidence=True,  # ignored
            show_cutoff=True,  # ignored
            sample_confidence=0.95,  # ignored
            colors=None,
            # ax=None,
        )
        plt.title(f"n_clusters{n_clusters}")
        plt.xlabel("lag time (steps)")  # 1 step = 1 interval in the trajectory
        plt.ylabel("implied timescale (steps)")  # 1 step = 1 interval in the trajectory
        plt.yscale("log")
        plt.savefig(
            f"{params.output_directory}/images/its_3d_n_clusters{n_clusters}_cut{params.cutoff}.png",
            bbox_inches="tight",
            pad_inches=0.1,
            dpi=300,
        )
        if params.show_picture:
            plt.show()
        plt.clf()
        plt.close()

        params.logger.info(f"Finished plotting ITS for n_clusters={n_clusters}")


# plot its
plot_its_3d(params=params)


# ## 3D: FEL along $d$
#
# - Objective:
#   Calculate the free energy landscape (FEL) along $d$ using all trials collectively.
#
# - Approach:
#   In contrast to the 1D analysis where each trial is processed separately, the 3D analysis is performed on the combined dataset of all trials.
#   Bayesian MSM generates a large number of samples, which allows us to estimate error bounds.
#   Only the mean FEL with its error bounds (mean $\pm$ error) is plotted, rather than all individual samples.
#
# 1. Determine the Bin Ranges:
#
#    - Let the cluster center of the $l$th cluster be denoted as $d^{l}$.
#    - Compute the overall minimum and maximum among all $d^{l}$ obtained from the clustering results.
#    - Divide the range between these minimum and maximum values into `params.n_bins` equal intervals.
#    - Define:
#      - The center of each interval as $d_i$.
#      - The length of each interval as $\delta d$.
#
# 2. Compute the FEL Using Bayesian MSM:
#
#    - Utilize the large number of samples generated by Bayesian MSM to estimate the FEL at each bin center $d_i$.
#    - For each $d_i$, compute the free energy as:
#      $$
#      \Delta W(d_i) = - RT \ln \left( \sum_{d_i - \frac{\delta d}{2} \le d^l \le d_i + \frac{\delta d}{2}} \pi_l \right)
#      $$
#      where $\pi_l$ is the stationary probability associated with the $l$ th cluster.
#    - From the ensemble of Bayesian MSM samples, derive the standard error of $\Delta W(d_i)$.
#
# 3. Plot the Averaged FEL with Error Bounds:
#
#    - Plot only the mean FEL with its error bounds, i.e., $\Delta W(d_i) \pm \text{error}$.
#    - Do not plot individual samples; only the aggregated statistics (mean $\pm$ error) are displayed.
def plot_fel_along_d_3d(
    params: Parameters,
    n_clusters: int,
    lag: int,
) -> None:
    # log
    params.logger.info(
        f"Starting plotting FEL along d for n_clusters={n_clusters} lag={lag}"
    )

    # load
    msm_result_csv = f"{params.output_directory}/result_csvs/3d_n_clusters{n_clusters}_lag{lag}_cut{params.cutoff}.csv"
    df = pd.read_csv(msm_result_csv)
    cluster_center_x = df["cluster_center_x"].values
    cluster_center_y = df["cluster_center_y"].values
    cluster_center_z = df["cluster_center_z"].values
    cluster_center_d = np.sqrt(
        cluster_center_x**2 + cluster_center_y**2 + cluster_center_z**2
    )
    stat_dists = df.iloc[:, 3:].values  # shape=(n_clusters, n_samples)

    # binning
    bins = np.linspace(
        start=min(cluster_center_d), stop=max(cluster_center_d), num=params.nbins
    )
    bins -= (bins[1] - bins[0]) / 2
    bin_centers = bins[:-1] + (bins[1] - bins[0]) / 2
    bin_energy_means = np.zeros(len(bins) - 1)
    bin_energy_stds = np.zeros(len(bins) - 1)
    for i in range(len(bins) - 1):
        bin_mask = (cluster_center_d >= bins[i]) & (cluster_center_d < bins[i + 1])
        pi_sum = np.sum(stat_dists[bin_mask, :], axis=0)
        energies = params._RT * np.log(pi_sum)
        bin_energy_means[i] = np.mean(energies)
        bin_energy_stds[i] = np.std(energies)
    bin_stat_dist_stes = bin_energy_stds / np.sqrt(len(stat_dists))
    bin_energy_means -= np.min(bin_energy_means)

    # plot
    plt.plot(
        bin_centers,
        bin_energy_means,
        color="tab:blue",
        linewidth=2,
    )
    plt.errorbar(
        bin_centers,
        bin_energy_means,
        yerr=bin_stat_dist_stes,
        fmt="o",
        label=r"Mean $\pm$SE",
        color="tab:blue",
        capsize=3,
    )
    plt.xlim(0, 13)
    plt.ylim(0, 20)
    plt.title(f"n_clusters={n_clusters} lag={lag}")
    plt.xlabel(r"$d$ [nm]")
    plt.ylabel(r"$\Delta W$ [kcal/mol]")
    plt.legend()
    plt.savefig(
        f"{params.output_directory}/images/fel_3d_n_clusters{n_clusters}_lag{lag}.png",
        bbox_inches="tight",
        pad_inches=0.1,
        dpi=300,
    )
    if params.show_picture:
        plt.show()
    plt.clf()
    plt.close()
    params.logger.info(
        f"Finished plotting FEL along d for n_clusters={n_clusters} lag={lag}"
    )


for n_clusters in n_clusters_3d_main:
    for lag in lags_3d_main:
        try:
            plot_fel_along_d_3d(params=params, n_clusters=n_clusters, lag=lag)
        except Exception:
            continue


def calc_vc_all_trials(
    params: Parameters,
    lower_unbound: float,
    upper_unbound: float,
) -> None:
    # check order of lower_unbound, upper_unbound
    if lower_unbound > upper_unbound:
        params.logger.error("lower_unbound must be less than upper_unbound. Exiting.")
        sys.exit(1)
    if upper_unbound > params.cutoff:
        params.logger.error("upper_unbound must be less than cut_off. Exiting.")
        sys.exit(1)

    # check if vc csv file exists
    vc_csv = f"{params.output_directory}/result_csvs/3d_VC.csv"
    if not Path(vc_csv).exists():
        with open(vc_csv, "w") as f:
            f.write("lower_unbound,upper_unbound,volume,correction\n")
    else:
        # check if already calculated
        df = pd.read_csv(vc_csv)
        if (
            len(
                df[
                    (df["lower_unbound"] == lower_unbound)
                    & (df["upper_unbound"] == upper_unbound)
                ]
            )
            > 0
        ):
            params.logger.debug(
                f"volume correction for lower_unbound={lower_unbound} upper_unbound={upper_unbound} already exists. Skipped"
            )
            return

    params.logger.info(
        f"Starting volume correction calculation for lower_unbound={lower_unbound} upper_unbound={upper_unbound}"
    )
    # calc vc
    inter_COM_vec_all = []
    for trial in params.n_trial_for_calc:
        a = np.load(f"{params.feature_3d_directory}/trial{trial}.npy")
        for each_rep in range(0, a.shape[0]):
            inter_COM_vec = a[each_rep]
            max_distance = np.max(np.linalg.norm(inter_COM_vec, axis=1))
            if max_distance > params.cutoff:
                continue
            inter_COM_vec_all.append(inter_COM_vec)
    inter_COM_vec_all = np.concatenate(inter_COM_vec_all)

    # select snapshots of unbound states
    inter_COM_d_all = np.linalg.norm(inter_COM_vec_all, axis=1)
    ix = np.where(
        (inter_COM_d_all >= lower_unbound) & (inter_COM_d_all < upper_unbound)
    )
    inter_COM_vec = inter_COM_vec_all[ix]

    # calculate volume correction
    if len(inter_COM_vec) > 0:
        hull = ConvexHull(inter_COM_vec)
        volume = hull.volume * 1000  # [angstrome^3]
        correction = params._RT * np.log(volume / 1661)  # [kcal/mol]
    else:
        volume = 0
        correction = 0

    # save to csv file
    with open(vc_csv, "a") as f:
        f.write(f"{lower_unbound},{upper_unbound},{volume},{correction}\n")

    params.logger.info(
        f"Finished volume correction calculation for lower_unbound={lower_unbound} upper_unbound={upper_unbound}"
    )


def calc_binding_energy_3d(
    params: Parameters,
    n_clusters: int,
    lag: int,
    lower_bound: float,
    upper_bound: float,
    lower_unbound: float,
    upper_unbound: float,
) -> None:
    # check the order of bound, unbound
    if lower_bound >= upper_bound:
        params.logger.error("lower_bound should be less than upper_bound. Exiting.")
        sys.exit(1)
    if lower_unbound >= upper_unbound:
        params.logger.error("lower_unbound should be less than upper_unbound. Exiting.")
        sys.exit(1)
    if upper_bound > lower_unbound:
        params.logger.error(
            "upper_bound should be less than or equal to lower_unbound. Exiting."
        )
        sys.exit(1)
    if upper_unbound > params.cutoff:
        params.logger.error(
            "upper_unbound should be less than or equal to cutoff. Exiting."
        )
        sys.exit(1)

    # check if csv file exists
    be_csv = f"{params.output_directory}/result_csvs/3d_binding_energy_summary.csv"
    if not Path(be_csv).exists():
        with open(be_csv, "w") as f:
            f.write(
                "n_clusters,lag,lower_bound,upper_bound,lower_unbound,upper_unbound,dG_PMF,dG_ste,VC\n"
            )
    else:
        # check if already calculated
        df = pd.read_csv(be_csv)
        if (
            len(
                df[
                    (df["n_clusters"] == n_clusters)
                    & (df["lag"] == lag)
                    & (df["lower_bound"] == lower_bound)
                    & (df["upper_bound"] == upper_bound)
                    & (df["lower_unbound"] == lower_unbound)
                    & (df["upper_unbound"] == upper_unbound)
                ]
            )
            > 0
        ):
            params.logger.debug(
                f"Binding energy for n_clusters={n_clusters} "
                f"lag={lag} "
                f"lower_bound={lower_bound} "
                f"upper_bound={upper_bound} "
                f"lower_unbound={lower_unbound} "
                f"upper_unbound={upper_unbound} "
                "already calculated. Skipping"
            )
            return

    # calc vc
    if params.do_volume_correction:
        vc_csv = f"{params.output_directory}/result_csvs/3d_VC.csv"
        calc_vc_all_trials(
            params=params,
            lower_unbound=lower_unbound,
            upper_unbound=upper_unbound,
        )

    # calc binding energy
    # load
    msm_result_csv = f"{params.output_directory}/result_csvs/3d_n_clusters{n_clusters}_lag{lag}_cut{params.cutoff}.csv"
    if not Path(msm_result_csv).exists():
        params.logger.info(
            f"MSM result csv file for n_clusters={n_clusters} lag={lag} doesn't exist. Skipping."
            f"First, run MSM."
        )
        return
    params.logger.info(
        f"Starting calculating binding energy for n_clusters={n_clusters} lag={lag} "
        f"lower_bound={lower_bound} upper_bound={upper_bound} lower_unbound={lower_unbound} upper_unbound={upper_unbound}"
    )
    df = pd.read_csv(msm_result_csv)
    cluster_center_x = df["cluster_center_x"].values
    cluster_center_y = df["cluster_center_y"].values
    cluster_center_z = df["cluster_center_z"].values
    cluster_center_d = np.sqrt(
        cluster_center_x**2 + cluster_center_y**2 + cluster_center_z**2
    )
    stat_dists = df.iloc[:, 3:].values  # shape=(n_clusters, n_samples)

    # select clusters in bound/unbound states
    ix_b = np.where(
        (cluster_center_d >= lower_bound) & (cluster_center_d < upper_bound)
    )[0]
    ix_u = np.where(
        (cluster_center_d >= lower_unbound) & (cluster_center_d < upper_unbound)
    )[0]
    pi_b = np.sum(stat_dists[ix_b, :], axis=0)
    pi_u = np.sum(stat_dists[ix_u, :], axis=0)
    dG = params._RT * np.log(pi_b / pi_u)  # [kcal/mol]
    dG_mean = np.mean(dG)
    dG_std = np.std(dG)
    dG_ste = dG_std / len(dG) ** 0.5

    # get vc information
    if params.do_volume_correction:
        vc_csv = f"{params.output_directory}/result_csvs/3d_VC.csv"
        vc_df = pd.read_csv(vc_csv)
        vc_df = vc_df[
            (vc_df["lower_unbound"] == lower_unbound)
            & (vc_df["upper_unbound"] == upper_unbound)
        ]
        vc = vc_df["correction"].values[0]
    else:
        vc = 0

    # save to csv file
    with open(be_csv, "a") as f:
        f.write(
            f"{n_clusters},{lag},{lower_bound},{upper_bound},{lower_unbound},{upper_unbound},{dG_mean},{dG_ste},{vc}\n"
        )

    params.logger.info(
        f"Finished calculating binding energy for n_clusters={n_clusters} lag={lag} "
        f"lower_bound={lower_bound} upper_bound={upper_bound} lower_unbound={lower_unbound} upper_unbound={upper_unbound}"
    )


for n_clusters in n_clusters_3d_main:
    for lag in lags_3d_main:
        calc_binding_energy_3d(
            params=params,
            n_clusters=n_clusters,
            lag=lag,
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            lower_unbound=lower_unbound,
            upper_unbound=params.cutoff,
        )


# ## 3D: $k_{on}, k_{off}$
# - calculate the rate constant $k_{on}, k_{off}$ using MFPTe
# - The relationship between MFPT and $k_{on}, k_{off}$ is as follows
#
# $$
# \begin{align}
# k_{off} &= \frac{1}{MFPT_{off}} \\
# k_{on} &= \frac{1}{MFPT_{on} C_{ligand}}
# \end{align}
# $$
#
# - where $MFPT_{off}$ is the mean first passage time(MFPT) from the bound state to the unbound state, $MFPT_{on}$ is MFPT from the unbound state to the bound state, and $C_{ligand}$ is the ligand concentration.
#
# - MFPT can be calculated from MSM theory. See [deeptime web site](https://deeptime-ml.github.io/latest/api/generated/deeptime.markov.tools.analysis.mfpt.html) for more information.
def calc_rate_constant_3d(
    params: Parameters,
    n_clusters: int,
    lag: int,
    lower_bound: float,
    upper_bound: float,
    lower_unbound: float,
    upper_unbound: float,
) -> None:
    # check if it already calculated
    rc_csv = f"{params.output_directory}/result_csvs/3d_rate_constant_summary.csv"
    if not Path(rc_csv).exists():
        with open(rc_csv, "w") as f:
            f.write(
                "n_clusters,lag,lower_bound,upper_bound,lower_unbound,upper_unbound,koff_mean,koff_ste,kon_mean,kon_ste\n"
            )
    else:
        # check if already calculated
        df = pd.read_csv(rc_csv)
        if (
            len(
                df[
                    (df["n_clusters"] == n_clusters)
                    & (df["lag"] == lag)
                    & (df["lower_bound"] == lower_bound)
                    & (df["upper_bound"] == upper_bound)
                    & (df["lower_unbound"] == lower_unbound)
                    & (df["upper_unbound"] == upper_unbound)
                ]
            )
            > 0
        ):
            params.logger.info(
                f"Rate constant for n_clusters={n_clusters} lag={lag} "
                f"lower_bound={lower_bound} upper_bound={upper_bound} "
                f"lower_unbound={lower_unbound} upper_unbound={upper_unbound} "
                f"already calculated. Skipping."
            )
            return

    # check the order of bound, unbound
    with open(
        f"{params.output_directory}/cluster_objs/3d_n_clusters{n_clusters}_cut{params.cutoff}.pkl",
        "rb",
    ) as f:
        cluster_pkl = pickle.load(f)
    with open(
        f"{params.output_directory}/MSM_objs/3d_n_clusters{n_clusters}_cut{params.cutoff}.pkl",
        "rb",
    ) as f:
        msm_pkl = pickle.load(f)
    with open(
        f"{params.output_directory}/Count_objs/3d_n_clusters{n_clusters}_cut{params.cutoff}.pkl",
        "rb",
    ) as f:
        count_pkl = pickle.load(f)

    if not (lag in count_pkl and lag in msm_pkl and lag in cluster_pkl):
        params.logger.info(f"lag={lag} does not exist. Skipping.")
        return

    cluster_centers = cluster_pkl["model"].cluster_centers
    largest_connected_set = deeptime.markov.tools.estimation.largest_connected_set(
        count_pkl[lag].count_matrix, directed=True
    )
    center_d = np.linalg.norm(cluster_centers, axis=1)[largest_connected_set]
    ix_u = np.where((center_d > lower_unbound) & (center_d < upper_unbound))[0]
    ix_b = np.where((center_d > lower_bound) & (center_d < upper_bound))[0]

    mfpt_on_list = []
    mfpt_off_list = []
    for i in range(len(msm_pkl[lag].samples)):
        T = msm_pkl[lag].samples[i].transition_matrix
        sd = msm_pkl[lag].samples[i].stationary_distribution
        mfpt_off = deeptime.markov.tools.analysis.mfpt(
            T=T, target=ix_u, origin=ix_b, tau=lag, mu=sd
        )
        mfpt_on = deeptime.markov.tools.analysis.mfpt(
            T=T, target=ix_b, origin=ix_u, tau=lag, mu=sd
        )
        mfpt_on_list.append(mfpt_on)
        mfpt_off_list.append(mfpt_off)

    mfpt_off_arr = np.array(mfpt_off_list)
    mfpt_on_arr = np.array(mfpt_on_list)

    koff_arr = 1 / mfpt_off_arr
    kon_arr = 1 / mfpt_on_arr / params.ligand_concentration

    koff_mean = np.mean(koff_arr)  # step^-1
    kon_mean = np.mean(kon_arr)  # step^-1 M^-1

    koff_ste = np.std(koff_arr) / len(koff_arr) ** 0.5
    kon_ste = np.std(kon_arr) / len(kon_arr) ** 0.5

    scaling_factor = 10**12 / params.dt  # [step] to [second]
    koff_mean *= scaling_factor  # s^-1
    kon_mean *= scaling_factor  # s^-1 M^-1
    koff_ste *= scaling_factor
    kon_ste *= scaling_factor

    # save
    with open(rc_csv, "a") as f:
        f.write(
            f"{n_clusters},{lag},{lower_bound},{upper_bound},{lower_unbound},{upper_unbound},{koff_mean},{koff_ste},{kon_mean},{kon_ste}\n"
        )

    # log
    params.logger.info(
        f"Finished calculating rate constants for n_clusters={n_clusters} lag={lag} "
        f"lower_bound={lower_bound} upper_bound={upper_bound} "
        f"lower_unbound={lower_unbound} upper_unbound={upper_unbound}"
    )


for n_clusters in n_clusters_3d_main:
    for lag in lags_3d_main:
        calc_rate_constant_3d(
            params=params,
            n_clusters=n_clusters,
            lag=lag,
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            lower_unbound=lower_unbound,
            upper_unbound=params.cutoff,
        )


# koff [s^-1]
# kon [s^-1 M^-1]


# ## 3D: FEL on 2D plane
# - plot the FEL on the x-y, x-z, and y-z planes.
def plot_fel_each_2d(
    params: Parameters,
    n_clusters: int,
    lag: int,
    coord_names: list[int],
) -> None:
    # load raw feature data
    features = []
    for trial in params.n_trial_for_calc:
        a = np.load(f"{params.feature_3d_directory}/trial{trial}.npy")
        for each_rep in range(0, a.shape[0]):
            features_per_rep = a[each_rep]
            max_distance = np.max(np.linalg.norm(features_per_rep, axis=1))
            if max_distance > params.cutoff:
                continue
            else:
                features.append(features_per_rep)
    features_concat = np.concatenate(features, axis=0)

    # load clustering obj
    clustering_model_pkl = f"{params.output_directory}/cluster_objs/3d_n_clusters{n_clusters}_cut{params.cutoff}.pkl"
    if not Path(clustering_model_pkl).exists():
        params.logger.info(
            f"Clustering obj for n_clusters={n_clusters} cutoff={params.cutoff} does not exist. Skipping. First, run clustering."
        )
        return
    with open(clustering_model_pkl, "rb") as f:
        clustering_result = pickle.load(f)
        dtrajs = clustering_result["trjs"]
        dtrajs_concat = np.concatenate(dtrajs, axis=0)

    # load msm
    msm_model_pkl = f"{params.output_directory}/MSM_objs/3d_n_clusters{n_clusters}_cut{params.cutoff}.pkl"
    if not Path(msm_model_pkl).exists():
        params.logger.info(
            f"MSM obj for n_clusters={n_clusters} cutoff={params.cutoff} does not exist. Skipping."
            f"First, run MSM."
        )
        return
    with open(msm_model_pkl, "rb") as f:
        msm_result = pickle.load(f)

    if lag not in msm_result:
        params.logger.info(f"lag={lag} does not exist. Skipping.")
        return

    # find weights for each snapshot in trajectory
    msm_model = msm_result[lag]
    weights = msm_model.prior.compute_trajectory_weights(dtrajs_concat)[0]

    # plot
    feature_dim = features[0].shape[-1]

    params.logger.info(
        f"Starting plotting 2D-FEL for n_clusters={n_clusters} lag={lag}"
    )
    for first_dim in range(feature_dim):
        for second_dim in range(first_dim + 1, feature_dim):
            first_coord = coord_names[first_dim]
            second_coord = coord_names[second_dim]

            energy_2d_obj = deeptime.util.energy2d(
                x=features_concat[:, first_dim],
                y=features_concat[:, second_dim],
                weights=weights,
                bins=100,
                kbt=-params._RT,  # negative sign to cancel the sign of energy
                shift_energy=True,
            )

            # plot
            plt.figure(figsize=(8, 6))
            ax, contour, cbar = deeptime.plots.plot_energy2d(
                energies=energy_2d_obj,
                ax=None,
                levels=100,
                contourf_kws=dict(cmap="nipy_spectral"),
                # contourf_kws=dict(cmap='nipy_spectral', vmin=0, vmax=20), # to limit the range of colorbar
                cbar=True,
                cbar_kws=None,
                cbar_ax=None,
            )
            cbar.set_label(r"$\Delta W$ [kcal/mol]")
            plt.title(
                f"n_clusters={n_clusters} lag={lag} dim={first_coord}-{second_coord}"
            )
            plt.xlabel(f"{first_coord} [nm]")
            plt.ylabel(f"{second_coord} [nm]")
            # plt.xlim(-10, 10) # to limit the range of x-axis
            # plt.ylim(-10, 10) # to limit the range of y-axis
            plt.savefig(
                f"{params.output_directory}/images/fel_2d_n_clusters{n_clusters}_lag{lag}_dim={first_coord}-{second_coord}_cut{params.cutoff}.png",
                bbox_inches="tight",
                pad_inches=0.1,
                dpi=300,
            )
            if params.show_picture:
                plt.show()
            plt.clf()
            plt.close()
    params.logger.info(
        f"Finished plotting 2D-FEL for n_clusters={n_clusters} lag={lag}"
    )


# run
for n_clusters in n_clusters_3d_main:
    for lag in lags_3d_main:
        plot_fel_each_2d(params, n_clusters, lag, coord_names=["X", "Y", "Z"])
