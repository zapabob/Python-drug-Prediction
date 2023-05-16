import numpy as np
from scipy.integrate import odeint
from scipy.integrate import trapz
import matplotlib.pyplot as plt

def model(C, t, k):
    dCdt = -k * C
    return dCdt

# Parameters
dose_mg = 80  # mg
BA = 0.98  # Bioavailability
dose = dose_mg * BA  # Adjusted for bioavailability
MW = 177.24  # molecular weight
dose_mol = dose / MW  # Convert dose to moles
bodyweight = 60  # kg
Vd = 36  # Distribution volume in L
dose_mol_per_L = dose_mol / Vd  # Initial concentration in mol/L
dose_mg_per_L = dose_mol_per_L * MW  # Initial concentration in mg/L

# Elimination rate constant
Cl = 1  # Clearance in L/h, assumed to be 1 L/h for this example
k = Cl / Vd  # per hour

# Time points
t = np.linspace(0, 6, 1000)  # 0 to 6 hours

# Solve ODE
C = odeint(model, dose_mg_per_L, t, args=(k,))

# Calculate AUC
AUC = trapz(C, t)

# Create figure
fig, ax = plt.subplots()

# Plot concentration over time
ax.plot(t, C, label='Concentration (mg/L)')
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Concentration (mg/L)')
ax.legend()

# Display AUC
print('AUC:', AUC)

# Save figure to desktop
fig.savefig('~/Desktop/pharmacokinetics.png')

# Show the plot
plt.show()
