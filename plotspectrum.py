import matplotlib.pyplot as plt
import numpy as np

def read_spectrum(filename):
    channels = []
    intensities = []
    with open(filename, "r") as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 2:
                try:
                    ch = int(parts[0])
                    intensity = float(parts[1])
                    channels.append(ch)
                    intensities.append(intensity)
                except ValueError:
                    continue
    return np.array(channels), np.array(intensities)

# 讀取光譜資料
spectrum_file = "spectrum_mpi.dat"
channels, intensities = read_spectrum(spectrum_file)

# 繪圖
plt.figure(figsize=(10, 5))
plt.plot(channels, intensities, color="blue", linewidth=1)
plt.xlabel("Channel")
plt.ylabel("Intensity")
plt.title("MPI-computed Spectrum")
plt.grid(True)
plt.tight_layout()
plt.show()
