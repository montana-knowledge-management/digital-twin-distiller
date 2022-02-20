for i in range(len(prod1)):
    plt.plot(range_c, (caser["torque"])[i], lw=2)
plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
plt.minorticks_on()
plt.xlabel("rotorangle [deg]")
plt.ylabel("Torque [Nm]")
plt.show()