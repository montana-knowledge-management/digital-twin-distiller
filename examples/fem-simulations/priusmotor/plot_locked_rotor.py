
import matplotlib.pyplot as plt
from digital_twin_distiller import setup_matplotlib
from digital_twin_distiller import ModelDir
import json

ModelDir.set_base(__file__)

x = []
y = []
alpha = []
T = []
with open(ModelDir.DATA / 'locked_rotor_50.json', 'r', encoding='utf-8') as f:
    res_ = json.load(f)
    alpha_50 = res_.pop('alpha')
    T_50 = res_.pop('T')

with open(ModelDir.DATA / 'locked_rotor_75.json', 'r', encoding='utf-8') as f:
    res_ = json.load(f)
    T_75 = res_.pop('T')

with open(ModelDir.DATA / 'locked_rotor_100.json', 'r', encoding='utf-8') as f:
    res_ = json.load(f)
    T_100 = res_.pop('T')

with open(ModelDir.DATA / 'locked_rotor_125.json', 'r', encoding='utf-8') as f:
    res_ = json.load(f)
    T_125 = res_.pop('T')

with open(ModelDir.DATA / 'locked_rotor_150.json', 'r', encoding='utf-8') as f:
    res_ = json.load(f)
    T_150 = res_.pop('T')

with open(ModelDir.DATA / 'locked_rotor_200.json', 'r', encoding='utf-8') as f:
    res_ = json.load(f)
    T_200 = res_.pop('T')

with open(ModelDir.DATA / 'locked_rotor_250.json', 'r', encoding='utf-8') as f:
    res_ = json.load(f)
    T_250 = res_.pop('T')

setup_matplotlib()

plt.figure(figsize=(12, 8))

plt.plot(alpha_50, T_50, lw=2)
plt.plot(alpha_50, T_75, lw=2)
plt.plot(alpha_50, T_100, lw=2)
plt.plot(alpha_50, T_125, lw=2)
plt.plot(alpha_50, T_150, lw=2)
plt.plot(alpha_50, T_200, lw=2)
plt.plot(alpha_50, T_250, lw=2)

plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
plt.minorticks_on()
plt.xlabel("Electrical angle [deg]")
plt.ylabel("Cogging Torque [Nm]")
# plt.legend()
# plt.savefig(ModelDir.MEDIA / "cogging_torque.pdf", bbox_inches="tight")
plt.savefig(ModelDir.MEDIA / "locked_rotor.png", bbox_inches="tight", dpi=650)
plt.show()

