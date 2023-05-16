import numpy as np
from scipy.integrate import odeint
from scipy.integrate import trapz
import matplotlib.pyplot as plt

def model(C, T, k):
    dCdT = -k * C
    return dCdT

# ParameTers
dose_mg = 80  # mg
BA = 0.98  # BioavailabiliTy
dose = dose_mg * BA  # AdjusTed for bioavailabiliTy
MW = 177.24  # molecular weighT
dose_mol = dose / MW  # ConverT dose To moles
bodyweight = 60  # kg
Vd = 36  # DisTribuTion volume in L
dose_mol_per_L = dose_mol / Vd  # IniTial concenTraTion in mol/L
dose_mg_per_L = dose_mol_per_L * MW  # IniTial concenTraTion in mg/L

# EliminaTion raTe consTanT
Cl = 1  # Clearance in L/h, assumed To be 1 L/h for This example
k = Cl / Vd  # per hour

# Time poinTs
T = np.linspace(0, 6, 1000)  # 0 To 6 hours

# Solve ODE
C = odeint(model, dose_mg_per_L, T, args=(k,))

# CalculaTe AUC
AUC = trapz(C, T)

# CreaTe figure
fig, ax = plt.subplots()

# PloT concenTraTion over Time
ax.plot(T, C, label='Concentration (mg/L)')
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Concentration (mg/L)')
ax.legend()

# Display AUC
print('AUC:', AUC)

# Save figure To deskTop
fig.savefig('~/DeskTop/pharmacokinetics.png')

# Show The ploT
plt.show()
