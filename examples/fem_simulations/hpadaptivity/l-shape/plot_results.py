import operator as op
from math import hypot
from pprint import pprint

import matplotlib.pyplot as plt
import pymongo
from matplotlib import colors
from model import gradu, u

from digital_twin_distiller import ModelDir, setup_matplotlib

ModelDir.set_base(__file__)
setup_matplotlib()

conn_str = "mongodb://localhost:27017"
client = pymongo.MongoClient(conn_str, serverSelectionTimeoutMS=5000, w=0)

db = client.lshape.results

COLOR_FEMM = "mediumblue"
COLOR_AGROS = "firebrick"
COLOR_COMSOL = "gold"
COLOR_NETGEN = "purple"

# V
def plot_v():
    plt.figure()

    r = tuple(db.aggregate([{"$match": {"solver": "femm"}}, {"$project": {"_id": 0, "V": "$V2.V", "elements": 1}}]))
    elements = tuple(map(op.itemgetter("elements"), r))
    v = tuple(map(op.itemgetter("V"), r))
    plt.plot(elements, v, "o-", color=COLOR_FEMM, markerfacecolor="white", markersize=3, label="FEMM")

    r = tuple(db.aggregate([{"$match": {"solver": "agros2d"}}, {"$project": {"_id": 0, "V": "$V2.V", "elements": 1}}]))
    elements = tuple(map(op.itemgetter("elements"), r))
    v = tuple(map(op.itemgetter("V"), r))
    plt.plot(elements, v, "o-", color=COLOR_AGROS, markerfacecolor="white", markersize=3, label="Agros2D")

    r = tuple(db.aggregate([{"$match": {"solver": "comsol35"}}, {"$project": {"_id": 0, "V": "$V2.V", "elements": 1}}]))
    elements = tuple(map(op.itemgetter("elements"), r))
    v = tuple(map(op.itemgetter("V"), r))
    plt.plot(elements, v, "o-", color=COLOR_COMSOL, markerfacecolor="white", markersize=3, label="Comsol 3.5")

    r = tuple(db.aggregate([{"$match": {"solver": "netgen"}}, {"$project": {"_id": 0, "V": "$V2.V", "elements": 1}}]))
    elements = tuple(map(op.itemgetter("elements"), r))
    v = tuple(map(op.itemgetter("V"), r))
    plt.plot(elements, v, "o-", color=COLOR_NETGEN, markerfacecolor="white", markersize=3, label="Ngsolve")

    ue = u(0.95, 0.95)
    plt.plot([10, 310000], [ue, ue], "k", label="exact")

    plt.xscale("log")
    plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xlabel("number of elements")
    plt.ylabel("Electric potential @ (0.95, 0.95) [V]")
    plt.legend()
    plt.savefig(ModelDir.MEDIA / "v.pdf", bbox_inches="tight")
    plt.savefig(ModelDir.MEDIA / "v.png", bbox_inches="tight", dpi=450)
    plt.show()


# ENERGY
def plot_energy():
    plt.figure()

    r = tuple(db.find({"solver": "femm"}, {"Energy": 1, "_id": 0, "elements": 1}))
    elements = tuple(map(op.itemgetter("elements"), r))
    energy = tuple(map(op.itemgetter("Energy"), r))
    plt.semilogx(elements, energy, "o-", color=COLOR_FEMM, markerfacecolor="white", markersize=3, label="FEMM")

    r = tuple(db.find({"solver": "agros2d"}, {"Energy": 1, "_id": 0, "elements": 1}))
    elements = tuple(map(op.itemgetter("elements"), r))
    energy = tuple(map(op.itemgetter("Energy"), r))
    plt.semilogx(elements, energy, "o-", color=COLOR_AGROS, markerfacecolor="white", markersize=3, label="Agros2D")

    plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xlabel("number of elements")
    plt.ylabel("Stored Energy [J]")
    plt.legend()
    plt.savefig(ModelDir.MEDIA / "energy.pdf", bbox_inches="tight")
    plt.savefig(ModelDir.MEDIA / "energy.png", bbox_inches="tight", dpi=450)
    plt.show()


# E0
def plot_e0():
    plt.figure()

    r = tuple(db.aggregate([{"$match": {"solver": "femm"}}, {"$project": {"_id": 0, "E": "$E0.abs", "elements": 1}}]))
    elements = tuple(map(op.itemgetter("elements"), r))
    e0 = tuple(map(op.itemgetter("E"), r))
    plt.semilogx(elements, e0, "o-", color=COLOR_FEMM, markerfacecolor="white", markersize=3, label="FEMM")

    r = tuple(
        db.aggregate([{"$match": {"solver": "agros2d"}}, {"$project": {"_id": 0, "E": "$E0.abs", "elements": 1}}])
    )
    elements = tuple(map(op.itemgetter("elements"), r))
    e0 = tuple(map(op.itemgetter("E"), r))
    plt.semilogx(elements, e0, "o-", color=COLOR_AGROS, markerfacecolor="white", markersize=3, label="Agros2D")

    plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xlabel("number of elements")
    plt.ylabel("|E| @ (0, 0)[V/m]")
    plt.legend()
    plt.savefig(ModelDir.MEDIA / "e0.pdf", bbox_inches="tight")
    plt.savefig(ModelDir.MEDIA / "e0.png", bbox_inches="tight", dpi=450)
    plt.show()


# E1
def plot_e1():
    plt.figure()

    r = tuple(db.aggregate([{"$match": {"solver": "femm"}}, {"$project": {"_id": 0, "E": "$E1.abs", "elements": 1}}]))
    elements = tuple(map(op.itemgetter("elements"), r))
    e1 = tuple(map(op.itemgetter("E"), r))
    plt.semilogx(elements, e1, "o-", color=COLOR_FEMM, markerfacecolor="white", markersize=3, label="FEMM")

    r = tuple(
        db.aggregate([{"$match": {"solver": "agros2d"}}, {"$project": {"_id": 0, "E": "$E1.abs", "elements": 1}}])
    )
    elements = tuple(map(op.itemgetter("elements"), r))
    e1 = tuple(map(op.itemgetter("E"), r))
    plt.semilogx(elements, e1, "o-", color=COLOR_AGROS, markerfacecolor="white", markersize=3, label="Agros2D")

    plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xlabel("number of elements")
    plt.ylabel("|E| @ (1e-9, 1e-9)[V/m]")
    plt.legend()
    plt.savefig(ModelDir.MEDIA / "e1.pdf", bbox_inches="tight")
    plt.savefig(ModelDir.MEDIA / "e1.png", bbox_inches="tight", dpi=450)
    plt.show()


if __name__ == "__main__":
    # print(gradu(0, 0))
    plot_v()
    # plot_energy()
    # plot_e0()
    # plot_e1()
