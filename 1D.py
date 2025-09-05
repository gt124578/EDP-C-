import numpy as np
import matplotlib.pyplot as plt

# --- Paramètres physiques et numériques ---
L = 1.0       # Longueur de la barre
T = 2.0       # Temps total
nx = 100      # Points en espace
D = 0.01      # Diffusivité
C = 0.03      # Vitesse d'advection
dx = L / (nx - 1)
dt = dx**2 / (2 * D)   # Condition CFL
nt = int(T / dt)

# Maillage spatial
x = np.linspace(0, L, nx)

# --- Condition initiale ---
u = np.zeros((nt, nx))
u[0, :] = np.sin(np.pi * x)

# --- Terme source ---
def function_f(x, tps):
    return (np.sin(np.pi * x)
            + D * (np.pi**2) * np.sin(np.pi * x) * (1 + tps)
            + C * np.pi * np.cos(np.pi * x) * (1 + tps))

# --- Conditions aux limites ---
def apply_boundary_conditions(u_n):
    u_n[0] = 0
    u_n[-1] = 0
    return u_n

# --- Boucle temporelle ---
for n in range(nt - 1):
    for i in range(1, nx - 1):
        convection = -C * (u[n, i + 1] - u[n, i - 1]) / (2 * dx)
        diffusion = D * (u[n, i + 1] - 2 * u[n, i] + u[n, i - 1]) / dx**2
        source = function_f(i * dx, n * dt)
        u[n + 1, i] = u[n, i] + dt * (convection + diffusion + source)
    u[n + 1, :] = apply_boundary_conditions(u[n + 1, :])


# --- Affichage du profil de température ---
plt.figure(figsize=(8,5))
for n in range(0, nt, nt // 5):  # 5 instants répartis
    plt.plot(x, u[n, :], label=f"t = {n*dt:.2f}")

plt.xlabel("Position x")
plt.ylabel("Température u(x,t)")
plt.title("Évolution de la température - Équation de la chaleur 1D")
plt.legend()
plt.grid(True)
plt.show()

# --- Carte spatio-temporelle ---
plt.figure(figsize=(8,5))
plt.imshow(u, extent=[0, L, T, 0], aspect='auto', cmap='hot')
plt.colorbar(label="Température")
plt.xlabel("Position x")
plt.ylabel("Temps t")
plt.title("Carte spatio-temporelle de la chaleur en 1D")
plt.show()