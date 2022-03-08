
import matplotlib.pyplot as plt
from digital_twin_distiller import setup_matplotlib
from digital_twin_distiller import ModelDir
import json
import pandas as pd
import numpy as np

ModelDir.set_base(__file__)

x = []
y = []
alpha = []
T = []
with open(ModelDir.DATA / 'locked_rotor_50.json', 'r', encoding='utf-8') as f:
    res_ = json.load(f)
    alpha_50 = res_.pop('alpha')
    T_50 = res_.pop('T')

with open(ModelDir.DATA / 'locked50.csv', 'r', encoding='utf-8') as f:
    res_ = pd.read_csv(f)
    alpha50 = res_["x"]
    T50 = res_["y"]

with open(ModelDir.DATA / 'locked_rotor_75.json', 'r', encoding='utf-8') as f:
    res_ = json.load(f)
    T_75 = res_.pop('T')

with open(ModelDir.DATA / 'locked75.csv', 'r', encoding='utf-8') as f:
    res_ = pd.read_csv(f)
    alpha75 = res_["x"]
    T75 = res_["y"]

with open(ModelDir.DATA / 'locked_rotor_100.json', 'r', encoding='utf-8') as f:
    res_ = json.load(f)
    T_100 = res_.pop('T')

with open(ModelDir.DATA / 'locked100.csv', 'r', encoding='utf-8') as f:
    res_ = pd.read_csv(f)
    alpha100 = res_["x"]
    T100 = res_["y"]

with open(ModelDir.DATA / 'locked_rotor_125.json', 'r', encoding='utf-8') as f:
    res_ = json.load(f)
    T_125 = res_.pop('T')

with open(ModelDir.DATA / 'locked125.csv', 'r', encoding='utf-8') as f:
    res_ = pd.read_csv(f)
    alpha125 = res_["x"]
    T125= res_["y"]

with open(ModelDir.DATA / 'locked_rotor_150.json', 'r', encoding='utf-8') as f:
    res_ = json.load(f)
    T_150 = res_.pop('T')

with open(ModelDir.DATA / 'locked150.csv', 'r', encoding='utf-8') as f:
    res_ = pd.read_csv(f)
    alpha150 = res_["x"]
    T150 = res_["y"]

with open(ModelDir.DATA / 'locked_rotor_200.json', 'r', encoding='utf-8') as f:
    res_ = json.load(f)
    T_200 = res_.pop('T')

with open(ModelDir.DATA / 'locked200.csv', 'r', encoding='utf-8') as f:
    res_ = pd.read_csv(f)
    alpha200 = res_["x"]
    T200 = res_["y"]

with open(ModelDir.DATA / 'locked_rotor_250.json', 'r', encoding='utf-8') as f:
    res_ = json.load(f)
    T_250 = res_.pop('T')

with open(ModelDir.DATA / 'locked250.csv', 'r', encoding='utf-8') as f:
    res_ = pd.read_csv(f)
    alpha250 = res_["x"]
    T250 = res_["y"]

setup_matplotlib()
a = 4
plt.figure(figsize=(6, 4))
z = np.polyfit(alpha_50, T_250, a)
predict = np.poly1d(z)
y = predict(alpha_50)
plt.plot(alpha_50, y, c="r", label="simulation (250A)")
z = np.polyfit(alpha50, T250, a)
predict = np.poly1d(z)
y = predict(alpha_50)
plt.plot(alpha_50, y, c="r", linestyle='--', label="measurement (250A)")
#plt.scatter(alpha_50, T_250, marker="o", c="g")
#plt.scatter(alpha50, T250, marker="x", c="r")

z = np.polyfit(alpha_50, T_200, a)
predict = np.poly1d(z)
y = predict(alpha_50)
plt.plot(alpha_50, y, c="g", label="simulation (200A)")
z = np.polyfit(alpha75, T200, a)
predict = np.poly1d(z)
y = predict(alpha_50)
plt.plot(alpha_50, y, c="g", linestyle='--', label="measurement (200A)")
#plt.scatter(alpha_50, T_200, lw=2, c="g")
#plt.scatter(alpha75, T200, lw=2, c="g")

z = np.polyfit(alpha_50, T_150, a)
predict = np.poly1d(z)
y = predict(alpha_50)
plt.plot(alpha_50, y, c="b", label="simulation (150A)")
z = np.polyfit(alpha75, T150, a)
predict = np.poly1d(z)
y = predict(alpha_50)
plt.plot(alpha_50, y, c="b", linestyle='--', label="measurement (150A)")
#plt.scatter(alpha_50, T_150, lw=2, c="b")
#plt.scatter(alpha75, T150, lw=2, c="b")

plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
plt.minorticks_on()
plt.xlabel("Electrical angle [deg]", fontsize=10)
plt.ylabel("Torque [Nm]", fontsize=10)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.legend(fontsize=10)
#plt.savefig(ModelDir.MEDIA / "cogging_torque.pdf", bbox_inches="tight")
plt.savefig(ModelDir.MEDIA / "locked_rotor.png", bbox_inches="tight", dpi=650)
#plt.show()

with open(ModelDir.DATA / 'cogging_torque.csv', 'r', encoding='utf-8') as f:
    res_ = pd.read_csv(f)
    x = res_["x"]
    y1 = res_["y1"]
    y2 = res_["y2"]
    y3 = res_["y3"]

plt.figure(figsize=(6, 4))
a = 4
z = np.polyfit(x, y1, a)
predict = np.poly1d(z)
y = predict(x)
plt.plot(x*4, y, c="r", label="case A")
z = np.polyfit(x, y2, a)
predict = np.poly1d(z)
y = predict(x)
plt.plot(x*4, y, c="g", linestyle='--', label="case B" )
z = np.polyfit(x, y3, a)
predict = np.poly1d(z)
y = predict(x)
plt.plot(x*4, y, c="b", linestyle='dotted', label="case C")


plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
plt.minorticks_on()
plt.xlabel("Electrical angle [deg]", fontsize=10)
plt.ylabel("Torque [Nm]", fontsize=10)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.legend(fontsize=10)
#plt.savefig(ModelDir.MEDIA / "cogging_torque.pdf", bbox_inches="tight")
plt.savefig(ModelDir.MEDIA / "cogging.png", bbox_inches="tight", dpi=650)
plt.show()