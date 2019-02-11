#!/usr/bin/env python3

"""Generate figures for instaGRAAL applied to Ectocarpus sp.

Load genome data and draw figures related to the reassembly of Ectocarpus
sp. (Baudry et al.)

"""

import numpy as np
from Bio import SeqIO
import matplotlib
from matplotlib import pyplot as plt

matplotlib.rc("font", serif="Helvetica Neue")
import seaborn as sns
import matplotlib.ticker as mtick
import pathlib

import scipy.stats

sns.set()

path_ecto_initial = pathlib.Path("data/genomes/ecto_sp_reference_genome.fa")
path_ecto_raw_instagraal = pathlib.Path(
    "instagraal_output/test_mcmc_4/genome.fasta"
)
path_ecto_polished = pathlib.Path("data/genomes/ecto_sp_polished_assembly.fa")

parameter_files = [
    pathlib.Path("instagraal_output/test_mcmc_4/" + u)
    for u in [
        "list_n_contigs.txt",
        "list_d_nuc.txt",
        "list_mean_len.txt",
        "list_likelihood.txt",
        "list_slope.txt",
        "list_fact.txt",
    ]
]

(N_CONTIGS, D_NUC, MEAN_LEN, LOG_LIKELIHOOD, SLOPE, PREFACTOR) = range(
    len(parameter_files)
)

lengths = sorted(
    (len(u.seq) for u in SeqIO.parse(path_ecto_raw_instagraal, "fasta")),
    reverse=True,
)

lengths_polished = sorted(
    (len(u.seq) for u in SeqIO.parse(path_ecto_polished, "fasta")),
    reverse=True,
)

(n_contigs, d_nuc, mean_len, likelihood, slope, fact) = parameter_datasets = [
    np.loadtxt(f) for f in parameter_files
]

label_parameters = [
    "Number of scaffolds",
    "Intra/inter distance threshold",
    "Mean length",
    "Log-likelihood",
    "Slope",
    "Pre-factor",
]


def plot_scaffold_sizes():

    plt.stem(
        range(1, 28),
        lengths[:27],
        markerfmt=" ",
        linefmt="red",
        label="Main scaffolds",
    )
    plt.stem(
        range(len(lengths) + 1)[28:],
        lengths[27:],
        markerfmt=" ",
        linefmt="blue",
        label="Remainder",
    )
    plt.yscale("log")
    plt.xlabel("Scaffold index")
    plt.ylabel("Scaffold size (base pairs)")
    plt.title("Raw instaGRAAL assembly scaffold sizes")

    plt.xlim(left=0.5, right=321.5)
    plt.hlines(
        [sum(lengths[28:])],
        xmin=27.5,
        xmax=321,
        colors="green",
        linestyles="dashed",
        label="Sum of the remainder",
    )
    handles, labels = plt.gca().get_legend_handles_labels()
    handles = [handles[1], handles[2], handles[0]]
    labels = [labels[1], labels[2], labels[0]]
    plt.legend(handles, labels)
    plt.show()
    plt.savefig(
        "raw_instagraal_assembly_scaffold_sizes.eps", bbox_inches="tight"
    )
    plt.close()

    plt.stem(
        range(1, 28),
        lengths_polished[:27],
        markerfmt=" ",
        linefmt="red",
        label="Main scaffolds",
    )
    plt.stem(
        range(len(lengths_polished) + 1)[28:],
        lengths_polished[27:],
        markerfmt=" ",
        linefmt="blue",
        label="Remainder",
    )
    plt.yscale("log")
    plt.xlabel("Scaffold index")
    plt.ylabel("Scaffold size (base pairs)")
    plt.title("Polished instaGRAAL assembly scaffold sizes")
    plt.xlim(left=0.5, right=len(lengths_polished) + 0.5)
    plt.hlines(
        [sum(lengths[28:])],
        xmin=27.5,
        xmax=len(lengths_polished),
        colors="green",
        linestyles="dashed",
        label="Sum of the remainder",
    )
    handles, labels = plt.gca().get_legend_handles_labels()
    handles = [handles[1], handles[2], handles[0]]
    labels = [labels[1], labels[2], labels[0]]
    plt.legend(handles, labels)
    plt.show()
    plt.savefig(
        "polished_instagraal_assembly_scaffold_sizes.eps", bbox_inches="tight"
    )
    plt.close()


def compute_iqr(dataset, length=10000):

    n = len(dataset)
    step = n // length

    iqrs = np.array(
        [scipy.stats.iqr(dataset[:i:step]) for i in range(0, n, step)],
        dtype=np.float64,
    )[1:]

    return iqrs


my_iqrs = [compute_iqr(dataset) for dataset in parameter_datasets]


def slugify(string):
    return string.replace("/", "_").replace(" ", "_").lower()


def generate_all_plots():
    for my_dataset, label in zip(parameter_datasets, label_parameters):
        sns.lineplot(data=my_dataset)
        plt.title("Evolution of the {}".format(label.lower()))
        plt.xlim(left=0, right=len(my_dataset) + 0.5)
        plt.xlabel("Iterations")
        plt.ylabel(label)
        plt.savefig("{}.png".format(slugify(label)), bbox_inches="tight")
        plt.savefig("{}.svg".format(slugify(label)), bbox_inches="tight")
        plt.close()

    for my_dataset, label in zip(my_iqrs, label_parameters):
        sns.lineplot(data=my_dataset)
        plt.title(
            "Evolution of the interquartile ranges of the {}".format(
                label.lower()
            )
        )
        plt.xlim(left=0, right=len(my_dataset) + 0.5)
        plt.xlabel("Cumulative iteration windows")
        plt.ylabel("IQR({})".format(label))
        plt.savefig("iqr_{}.png".format(slugify(label)), bbox_inches="tight")
        plt.savefig("iqr_{}.svg".format(slugify(label)), bbox_inches="tight")
        plt.close()


def generate_size_dist():
    fig = plt.figure()
    ax3 = fig.add_subplot()
    ax3.stem(
        range(1, 28),
        lengths[:27],
        markerfmt=" ",
        linefmt="red",
        label="Main scaffolds",
    )
    ax3.stem(
        range(len(lengths) + 1)[28:],
        lengths[27:],
        markerfmt="",
        linefmt="blue",
        label="Remainder",
    )
    ax3.hlines(
        [sum(lengths[28:])],
        xmin=27.5,
        xmax=321,
        colors="green",
        linestyles="dashed",
        label="Sum of the remainder",
    )
    handles, labels = ax3.get_legend_handles_labels()
    handles = [handles[1], handles[2], handles[0]]
    labels = [labels[1], labels[2], labels[0]]
    ax3.legend(handles, labels)
    ax3.set_yscale("log")
    ax3.set_xlabel("Scaffold index")
    ax3.set_ylabel("Scaffold size (bp)")
    ax3.set_title("InstaGRAAL assembly scaffold sizes")
    plt.show()


def plot_matrices_and_parameters():
    fig = plt.figure(figsize=(4, 10))
    ax4 = fig.add_subplot(4, 5, 5)
    ax5 = fig.add_subplot(4, 5, 10, sharex=ax4)
    ax6 = fig.add_subplot(4, 5, 15, sharex=ax4)
    ax7 = fig.add_subplot(4, 5, 20, sharex=ax4)

    ax4.plot(parameter_datasets[4], color="red")
    ax4.set_ylabel("Exponent")
    plt.setp(ax4.get_xticklabels(), visible=False)
    ax5.plot(parameter_datasets[5], color="red")
    ax5.set_ylabel("Pre-factor")
    plt.setp(ax5.get_xticklabels(), visible=False)
    ax6.plot(parameter_datasets[1], color="red")
    ax6.set_ylabel("Mean trans contacts")
    plt.setp(ax6.get_xticklabels(), visible=False)
    ax7.plot(parameter_datasets[3], color="red")
    ax7.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.0e"))
    ax7.set_ylabel("Log-likelihood")
    ax7.set_xlabel("Iterations")

    for ax in ax4, ax5, ax6, ax7:
        ax.set_xticks(np.round(np.linspace(0, 800000, 2), 2))

    fig.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)
    plt.show()


matplotlib.rc("font", family="sans-serif")
matplotlib.rc("font", serif="Helvetica")
matplotlib.rc("text", usetex="false")
matplotlib.rcParams.update({"font.size": 14})

plot_scaffold_sizes()
generate_size_dist()
plot_matrices_and_parameters()
