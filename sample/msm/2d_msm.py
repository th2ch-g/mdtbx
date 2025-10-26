import deeptime
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
import re
import pickle

n_clusters_list = [50, 100, 200, 300, 400, 500]
lags_list = [1, 10, 100, 1000, 5000]
T = 310
RT = 8.314 * T / 1000 / 4.184

features = []
open_tmscores = []
close_tmscores = []
for i in Path(".").glob("*_close.tmscore"):
    tmp_close = []
    with open(str(i)) as ref:
        for idx, line in enumerate(ref):
            line = line.rstrip()
            if idx == 0:
                continue
            a = re.findall(r"\S+", line)
            tmp_close.append(float(a[3]))

    j = str(i).replace("close", "open")
    tmp_open = []
    with open(j) as ref:
        for idx, line in enumerate(ref):
            line = line.rstrip()
            if idx == 0:
                continue
            a = re.findall(r"\S+", line)
            tmp_open.append(float(a[3]))

    tmp = np.array([tmp_open, tmp_close])
    tmp = tmp.T
    open_tmscores.extend(tmp_open)
    close_tmscores.extend(tmp_close)

    features.append(tmp)

open_tmscores = np.array(open_tmscores)
close_tmscores = np.array(close_tmscores)


print(np.array(features).shape)


if Path("msm_result.pkl").is_file():
    with open("msm_result.pkl", "rb") as f:
        msm_result = pickle.load(f)
else:
    msm_result = dict()

if Path("clustering_trjs_result.pkl").is_file():
    with open("clustering_trjs_result.pkl", "rb") as f:
        clustering_trjs_result = pickle.load(f)
else:
    clustering_trjs_result = dict()

if Path("clustering_model_result.pkl").is_file():
    with open("clustering_model_result.pkl", "rb") as f:
        clustering_model_result = pickle.load(f)
else:
    clustering_model_result = dict()

if Path("count_result.pkl").is_file():
    with open("count_result.pkl", "rb") as f:
        count_result = pickle.load(f)
else:
    count_result = dict()

# clustering
for n_clusters in n_clusters_list:
    if n_clusters in clustering_trjs_result.keys():
        print(f"skip {n_clusters}")
        clustering_trjs = clustering_trjs_result[n_clusters]
        clustering_model = clustering_model_result[n_clusters]
    else:
        print(f"start clustering {n_clusters}")
        estimator = deeptime.clustering.KMeans(
            n_clusters=n_clusters,
            max_iter=2000,
            metric="euclidean",
            tolerance=1e-05,
            init_strategy="kmeans++",
            fixed_seed=False,
            n_jobs=None,
            initial_centers=None,
            progress=None,
        )

        print("estimator fit")
        estimator.fit(np.concatenate(features))

        clustering_trjs = []  # trjs
        for feature in features:
            clustering_trjs.append(estimator.transform(feature))
        print("estimator fetch model")
        clustering_model = estimator.fetch_model()  # models

        clustering_trjs_result[n_clusters] = clustering_trjs
        clustering_model_result[n_clusters] = clustering_model

    count_result_lag = dict()
    msm_result_lag = dict()
    for lag in lags_list:
        if n_clusters in msm_result.keys() and lag in count_result[n_clusters].keys():
            print(f"skip {lag}")
            msm_result_lag[lag] = msm_result[n_clusters][lag]
            count_result_lag[lag] = count_result[n_clusters][lag]
        else:
            print(f"start lag {lag}")
            count_estimator = deeptime.markov.TransitionCountEstimator(
                lagtime=lag,
                count_mode="effective",
                n_states=None,
                sparse=False,
            )

            print("count estimator fit")
            count_model = count_estimator.fit(
                clustering_trjs, n_jobs=None
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
                print("msm estimator fit")
                msm_model = msm_estimator.fit(count_model).fetch_model()
            except Exception as e:
                print("Error occured")
                print(e)
                continue

            # save
            msm_result_lag[lag] = msm_model
            count_result_lag[lag] = count_model

            print("end lag")

    msm_result[n_clusters] = msm_result_lag
    count_result[n_clusters] = count_result_lag

    # plot its
    msm_models = []
    for lag in sorted(msm_result[n_clusters].keys()):
        msm_models.append(msm_result[n_clusters][lag])

    its_data = deeptime.util.validation.implied_timescales(
        models=msm_models,
        n_its=10,
    )

    plt.figure(figsize=(8, 6))
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
    plt.tight_layout()
    plt.savefig(
        f"its_2d_n_clusters{n_clusters}.png",
        bbox_inches="tight",
        pad_inches=0.1,
        dpi=300,
        transparent=True,
    )
    plt.savefig(
        f"its_2d_n_clusters{n_clusters}.pdf",
        bbox_inches="tight",
        pad_inches=0.1,
        dpi=300,
        transparent=True,
    )
    plt.clf()
    plt.close()


with open("msm_result.pkl", "wb") as f:
    pickle.dump(msm_result, f)

with open("clustering_trjs_result.pkl", "wb") as f:
    pickle.dump(clustering_trjs_result, f)

with open("clustering_model_result.pkl", "wb") as f:
    pickle.dump(clustering_model_result, f)

with open("count_result.pkl", "wb") as f:
    pickle.dump(count_result, f)


for n_clusters in n_clusters_list:
    print(f"n_clusters: {n_clusters}")
    for lag in lags_list:
        print(f"lag: {lag}")
        msm_model = msm_result[n_clusters][lag]
        clustering_trjs = clustering_trjs_result[n_clusters]
        dtrajs_concat = np.concatenate(clustering_trjs, axis=0)
        weights = msm_model.prior.compute_trajectory_weights(dtrajs_concat)[0]

        plt.figure(figsize=(8, 6))
        energy2d = deeptime.util.energy2d(
            x=open_tmscores,
            y=close_tmscores,
            kbt=RT,
            weights=weights,
            # weights=None, # None means same weight for all samples
            shift_energy=True,
        )
        ax, contour, cbar = deeptime.plots.plot_energy2d(
            energy2d,
            ax=None,
            levels=100,
            cbar=True,
            cbar_kws=None,
            cbar_ax=None,
        )
        # plt.scatter(open_tmscores, close_tmscores, alpha=0.1, marker=".", s=1)
        cbar.set_label(r"$\Delta W$ [kcal/mol]")
        plt.grid(True)
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.xlabel("Similarity to open structures [TM-score]")
        plt.ylabel("Similarity to close structures [TM-score]")
        plt.tight_layout()
        plt.savefig(
            f"ene2d_n_clusters{n_clusters}_lag{lag}.png", dpi=300, transparent=True
        )
        plt.savefig(
            f"ene2d_n_clusters{n_clusters}_lag{lag}.pdf", dpi=300, transparent=True
        )
        plt.clf()
        plt.close()

        print("end clustering")
