import numpy as np
import matplotlib.pyplot as plt

# Range of u in (0, 1]
u = np.linspace(0, 1, 200)[1:]

# Some d_F values
d_values = [10, 8, 6, 4, 2, 0]

plt.figure(figsize=(6, 4), dpi=300)

# Calculate the alpha values
alphas = np.linspace(0.4, 1, len(d_values))

for d, alpha in zip(d_values, alphas):
    # Since D_{1,G} = D_{2,G} = 1:
    # eta_G(d_F, u) = d_F*(u + sqrt(1-u^2)) when u>0
    # we multiply by u => eta_G(d_F, u)*u = d_F*(u + sqrt(1-u^2))*u
    y = -d * (u + np.sqrt(1 - u**2)) * u
    plt.plot(u, y, label=f"d_F = {d}", color='#65107a', alpha=alpha)

plt.title(r"$\partial_t d_F = -\eta_G(d_F, u)\,u = -d_F\,u\,(u + \sqrt{1 - u^2})$")
plt.xlabel(r"$u$")
plt.ylabel(r"$\partial_t d_F$")
plt.xlim(left=-0.2, right=1.1)
tick_line_width = 1
plt.axvline(x=1, color='gray', linestyle='--', alpha=0.7, linewidth=tick_line_width)
plt.axhline(y=0, color='black', alpha=0.7, linewidth=tick_line_width)
plt.axvline(x=0, color='black', alpha=0.7, linewidth=tick_line_width)

# Calculate and plot the minima for each d_value
for d in d_values:
    y_min = -d * (1 + np.sqrt(2)) / 2
    plt.axhline(y=y_min, color='gray', linestyle='--', alpha=0.7, linewidth=tick_line_width)
    plt.text(-0.05, y_min, f"{y_min:.2f}", va='center', ha='right', color='gray', alpha=0.7, bbox=dict(facecolor='white', edgecolor='none', pad=2))

plt.legend(loc='lower left', bbox_to_anchor=(0.17, 0.01))
plt.xticks([])
plt.yticks([])
plt.show()
