import matplotlib.pyplot as plt
import sys

X = "10"

if len(sys.argv) < 2:
    print("python plotratio.py X")
    sys.exit(1)

X = sys.argv[1]

def read_data(filename):
    times = []
    with open(filename, "r") as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) >= 1:
                try:
                    times.append(float(parts[0]))
                except ValueError:
                    continue
    return times



result1 = read_data("result1.dat")
resultX = read_data("result" + X + ".dat")
#resultX_subset = read_data("result" + X + "_subset.dat")


x1 = list(range(1, len(result1) + 1))
xx = list(range(1, len(resultX) + 1))
#xs = list(range(1, len(resultX_subset) + 1))

# avg result1
avg_result1 = sum(result1) / len(result1)
# Plot ratio between resultX and avg_result1
resultX_ratio = [avg_result1 / x for x in resultX]
#resultX_subset_ratio = [avg_result1 / x for x in resultX_subset]

plt.figure(figsize=(8, 5))
plt.plot(xx, resultX_ratio, label="T(non_reentrant) / T(fits_read_img)", marker='o', color='orange')
#plt.plot(xs, resultX_subset_ratio, label="T(non_reentrant) / T(fits_read_subset)", marker='s', color='blue')

plt.xlabel("Trial")
plt.ylabel("Speed Ratio (non_reentrant: 1)")
plt.title("FITS Read Time Comparison")
plt.ylim(top=11)
plt.ylim(bottom=0)
plt.legend()
plt.grid()

plt.savefig("1v"+X+"_mpi_alma16G_1b1_ratio_g++-15_c++20_lmpi.png", dpi=300)
plt.show()
