import matplotlib.pyplot as plt


agrosdata = {"F1": [], "F2": [], "F3": [], "nodes": [], "r": []}

femmdata = {"F1": [], "F2": [], "F3": [], "nodes": [], "r": []}

with open("f_conv_stats.csv") as f:
    for line in f:
        record = line.strip().split(",")
        key = record.pop(0)
        record = [float(ri) for ri in record]
        if key == "agros2d":
            agrosdata["F1"].append(record.pop(0))
            agrosdata["F2"].append(record.pop(0))
            agrosdata["F3"].append(record.pop(0))
            agrosdata["nodes"].append(record.pop(0))
            agrosdata["r"].append(record.copy())

        if key == "femm":
            femmdata["F1"].append(record.pop(0))
            femmdata["F2"].append(record.pop(0))
            femmdata["F3"].append(record.pop(0))
            femmdata["nodes"].append(record.pop(0))
            femmdata["r"].append(record.copy())

# plt.figure()
# plt.scatter(agrosdata["F1"], agrosdata["F2"])
# plt.show()
plt.figure()
plt.scatter(agrosdata["nodes"], agrosdata["F1"], label="Agros")
plt.scatter(femmdata["nodes"], femmdata["F1"], label="Femm")
plt.legend()
plt.show()

print(len(agrosdata["nodes"]))
print(len(femmdata["nodes"]))
