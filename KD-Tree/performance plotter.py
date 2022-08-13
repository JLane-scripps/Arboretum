import ast
import os.path
from matplotlib import pyplot as plt

path = 'Sim Log'
folder = os.scandir(path)
labels = []
y_add = []
y_search = []

for file in folder:
    x = []
    if file.is_file():
        with open(file, 'r') as log:
            f_name = log.name.split('\\')
            labels.append(f_name[1])
            performance_dict = ast.literal_eval(log.read())
            plt.plot(list(performance_dict.keys()), [performance_dict[key]['add_time'] for key in performance_dict],
                     label=f_name[1].split(".")[0])

plt.xlabel("Num PSMs")
plt.ylabel("Time")
plt.title("PSM Tree Add Time")
plt.legend()
plt.show()

folder = os.scandir(path)
for file in folder:
    x = []
    if file.is_file():
        with open(file, 'r') as log:
            f_name = log.name.split('\\')
            labels.append(f_name[1])
            performance_dict = ast.literal_eval(log.read())
            plt.plot(list(performance_dict.keys()), [performance_dict[key]['search_time'] for key in performance_dict],
                     label=f_name[1].split(".")[0])

plt.xlabel("Num PSMs")
plt.ylabel("Time")
plt.title("PSM Tree Search Time")
plt.legend()
plt.show()

print("add times are", y_add)
print("search times are", y_search)
print(x)
print(labels)
folder.close()
file_name = os.path.join('Arboretum', "point_data.txt")


