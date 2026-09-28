import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# 1. Load the TRANSP NetCDF file
file_path = "174082C01.CDF"
nc = Dataset(file_path, "r")


# 2. Extract coordinates (Time and Radial Grid)
# TIME is usually 1D. TIME_V denotes the time axis for profiles.
time = nc.variables["TIME"][:]  # Time array in seconds
rhotor = nc.variables["X"][:] #rhotor ?
dvol = nc.variables["DVOL"][:]  # Volume profile in cm^3

# 4. Pick a specific time slice to plot
# Let's find the index closest to a target time (e.g., 2.5 seconds)
#target_time = 2.5  
#time_idx = np.abs(time - target_time).argmin()
time_idx = -1 
actual_time = time[time_idx]

# Slice data for this specific time
rhotor_slice = rhotor[time_idx, :] if rhotor.ndim == 2 else rhotor
dvol_slice = dvol[time_idx,:] 

# 3. Extract Transport Powers (Conduction + Convection)
# Units in TRANSP NetCDF are typically Watts. We will convert to MW (/1e6).
# Electron channel
p_cond_e = nc.variables["PCNDE"][time_idx,:]  # Electron conduction power
p_conv_e = nc.variables["PCNVE"][time_idx,:]  # Electron convection power
p_total_e = (p_cond_e + p_conv_e)

# Ion channel
p_cond_i = nc.variables["PCOND"][time_idx,:]  # Ion conduction power
p_conv_i = nc.variables["PCONV"][time_idx,:]  # Ion convection power
p_total_i = (p_cond_i + p_conv_i)

# Perform the Volume Integration to get Watts, then divide by 1e6 for MW
p_total_e_mw = np.cumsum(p_total_e*dvol_slice) / 1e6
p_total_i_mw = np.cumsum(p_total_i*dvol_slice) / 1e6




# 5. Plot the Results
plt.figure(figsize=(8, 5))
plt.plot(rhotor_slice, p_total_e_mw, label="Electron Channel ($P_{cond} + P_{conv}$)", color="crimson", lw=2)
plt.plot(rhotor_slice, p_total_i_mw, label="Ion Channel ($P_{cond} + P_{conv}$)", color="royalblue", lw=2)

# Add explicit markers at x = 0.85
# 8. Find the values at the coordinate index closest to 0.85
target_x = 0.85
idx_085 = np.abs(rhotor_slice - target_x).argmin()
actual_x = rhotor_slice[idx_085]
val_e = p_total_e_mw[idx_085]
val_i = p_total_i_mw[idx_085]

plt.scatter([actual_x], [val_e], color="crimson", s=40, zorder=5)
plt.scatter([actual_x], [val_i], color="royalblue", s=40, zorder=5)

# Annotate values with slight offsets to ensure they remain clear and legible
plt.annotate(f"{val_e:.2f} MW", (actual_x, val_e), textcoords="offset points", xytext=(-45, 15),
             ha='center', color="crimson", weight='bold',
             arrowprops=dict(arrowstyle="->", color="crimson", lw=1))

plt.annotate(f"{val_i:.2f} MW", (actual_x, val_i), textcoords="offset points", xytext=(-45, -20),
             ha='center', color="royalblue", weight='bold',
             arrowprops=dict(arrowstyle="->", color="royalblue", lw=1))

# Draw a vertical reference line at the target coordinate
plt.axvline(x=actual_x, color="gray", linestyle=":", alpha=0.7, label=f"X = {actual_x:.2f}")

plt.title(f"TRANSP Run 174082C01 - Transport Power Channels at t = {actual_time:.3f} s")
plt.xlabel(r"Normalized Radius ($\rho_{tor}$)")
plt.ylabel("Power (MW)")
plt.xlim(0, 1)
plt.grid(True, linestyle="--", alpha=0.5)
plt.legend()

# Display the plot
plt.tight_layout()
plt.show()
