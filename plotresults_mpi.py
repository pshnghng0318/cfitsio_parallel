import matplotlib.pyplot as plt
import sys

X = "10"
if len(sys.argv) < 2:
    print("python plotresults_mpi.py X")
    sys.exit(1)

X = sys.argv[1]
print("X =", X)

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

x1 = list(range(1, len(result1) + 1))
xx = list(range(1, len(resultX) + 1))

plt.figure(figsize=(8, 5))
plt.plot(x1, result1, marker="o", linestyle="-", label="Single Thread (result1)")
plt.plot(xx, resultX, marker="s", linestyle="--", label="MPI " + X + " CPU (result" + X + ")")

plt.xlabel("Trial")
plt.ylabel("Time (sec)")
plt.title("FITS Read Time Comparison")
plt.ylim(top=60)
plt.ylim(bottom=0)
plt.legend()
plt.grid()

plt.savefig("1v"+X+"_mpi_alma16G_1b1_g++-15_c++20_lmpi.png", dpi=300)
plt.show()
