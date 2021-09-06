import csv
import matplotlib.pyplot as plt
from pathlib import Path


def setup_matplotlib():
    plt.style.use(["default", "seaborn-bright"])
    mm2inch = lambda x: 0.03937007874 * x
    plt.rcParams["figure.figsize"] = mm2inch(140), mm2inch(78.75)
    plt.rcParams["lines.linewidth"] = 1

    SMALL_SIZE = 8
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 14

    plt.rc("font", size=MEDIUM_SIZE)  # controls default text sizes
    plt.rc("axes", titlesize=MEDIUM_SIZE)  # fontsize of the axes title
    plt.rc("axes", labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc("xtick", labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
    plt.rc("ytick", labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
    plt.rc("legend", fontsize=MEDIUM_SIZE)  # legend fontsize
    plt.rc("figure", titlesize=MEDIUM_SIZE)  # fontsize of the figure title

def plot_torque():
    setup_matplotlib()
    dir_base = Path(__file__).parent
    dir_data = dir_base / "data"
    dir_media = dir_base / "media"
    dir_ref = dir_data / "reference_torque"
    dir_sim = dir_data / "torque"

    getlabel = lambda f: int(f.stem.split("_").pop())

    files_ref = list(dir_ref.glob("torque_*.csv"))
    files_ref.sort(key=getlabel)

    files_sim = list(dir_sim.glob("torque_*.csv"))
    files_sim.sort(key=getlabel)
    T_max_sims = []
    T_max_refs = []
    I = []
    for file_ref, file_sim in zip(files_ref, files_sim):
        label_ref = label = getlabel(file_ref)
        label_sim = label = getlabel(file_sim)
        assert label_sim==label_ref

        theta_sim = []
        T_sim = []
        theta_ref = []
        T_ref = []

        with open(file_ref, 'r') as f_ref, open(file_sim, 'r') as f_sim:
            for theta_i, Ti in csv.reader(f_ref, quoting=csv.QUOTE_NONNUMERIC):
                theta_ref.append((theta_i-90)*4)
                T_ref.append(Ti)

            for theta_i, Ti in csv.reader(f_sim, quoting=csv.QUOTE_NONNUMERIC):
                theta_sim.append(theta_i)
                T_sim.append(Ti*8)

        max_T_ref = max(T_ref)
        max_T_sim = max(T_sim)
        theta_max_ref = theta_ref[T_ref.index(max_T_ref)]
        theta_max_sim = theta_sim[T_sim.index(max_T_sim)]

        I.append(label)
        T_max_sims.append(max_T_sim)
        T_max_refs.append(max_T_ref)

        p = plt.plot(theta_sim, T_sim, label=f"{label} A")
        ci = p[-1].get_color()
        plt.plot(theta_ref, T_ref, color=ci, linestyle='--')

        plt.scatter(theta_max_ref, max_T_ref, c=ci)
        plt.scatter(theta_max_sim, max_T_sim, c=ci)

        print(f"I = {label:>4} A\t {max_T_ref=:>6.2f} Nm\t {max_T_sim=:>6.2f} Nm\t"
              f"{theta_max_ref=:.1f} °\t {theta_max_sim=:.1f} °")


    plt.grid(b=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(b=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xlabel("Electrical angle [°]")
    plt.ylabel("Torque [Nm]")
    plt.legend()
    plt.savefig(dir_media / "torque.png", dpi=550, bbox_inches="tight")
    plt.show()


    plt.figure()
    plt.plot(I, T_max_sims, '-o', label='Simulation')
    plt.plot(I, T_max_refs, '-o', label='Reference')
    plt.grid(b=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(b=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xlabel("Current [A]")
    plt.ylabel("Maximum Torque [Nm]")
    plt.legend()
    plt.savefig(dir_media / "current_max_torque.png", dpi=500, bbox_inches="tight")
    plt.show()




if __name__=='__main__':
    # TODO: ref interpolation
    plot_torque()