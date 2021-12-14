
import matplotlib.pyplot as plt
from digital_twin_distiller import setup_matplotlib
from digital_twin_distiller import ModelDir
import json

ModelDir.set_base(__file__)

alpha = []
T = []
with open(ModelDir.DATA / 'locked_rotor.json', 'r', encoding='utf-8') as f:
    res_ = json.load(f)
    alpha = res_.pop('alpha')
    T = res_.pop('T')


setup_matplotlib()

plt.figure(figsize=(12, 8))

plt.plot(alpha, T, lw=2)
plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
plt.minorticks_on()
plt.xlabel("Electrical angle [deg]")
plt.ylabel("Cogging Torque [Nm]")
# plt.legend()
# plt.savefig(ModelDir.MEDIA / "cogging_torque.pdf", bbox_inches="tight")
plt.savefig(ModelDir.MEDIA / "locked_rotor.png", bbox_inches="tight", dpi=650)
plt.show()

