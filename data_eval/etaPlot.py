import numpy as np
import matplotlib.pyplot as plt

# Range of u in [-1, 1]
u = np.linspace(-1, 1, 400)

# d_F value
d_F = 5

# Calculate eta_G(d_F, u)
y = d_F * (np.abs(u) + np.sqrt(np.clip(1 - u**2, 0, None)))

xMax1 = -np.sqrt(2) / 2
xMax2 = np.sqrt(2) / 2
yMax = d_F * (np.abs(xMax1) + np.sqrt(np.clip(1 - xMax1**2, 0, None)))

# Clamp yMax to the desired range
yMax = np.clip(yMax, -5, 9)

plt.figure(figsize=(6, 4), dpi=300)
plt.plot(u, y, color='#65107a')
# Calculate eta_G(d_F, u) for u < 0
y_dashed = d_F * (u + np.sqrt(np.clip(1 - u**2, 0, None)))

# Plot the dashed line for u < 0
plt.plot(u[u < 0], y_dashed[u < 0], linestyle='--', color='#65107a', alpha=0.6)

plt.title(r"$\eta_G(d_F, u) = d_F\,(\left|u\right| + \sqrt{1 - u^2})$")
plt.xlabel(r"$u$")
plt.ylabel(r"$\partial_t G$")
plt.ylim(bottom=-5, top=9)

# Mark the specified values with dashed grid lines
tick_line_width = 1
plt.axvline(x=-1, color='gray', linestyle='--', alpha=0.7, linewidth=tick_line_width)
plt.axvline(x=1, color='gray', linestyle='--', alpha=0.7, linewidth=tick_line_width)
plt.axvline(x=xMax1, color='gray', linestyle='--', alpha=0.7, linewidth=tick_line_width)
plt.axvline(x=xMax2, color='gray', linestyle='--', alpha=0.7, linewidth=tick_line_width)
plt.axhline(y=d_F, color='gray', linestyle='--', alpha=0.7, linewidth=tick_line_width)
plt.axhline(y=yMax, color='gray', linestyle='--', alpha=0.7, linewidth=tick_line_width)
plt.axhline(y=0, color='black', alpha=0.7, linewidth=tick_line_width)
plt.axvline(x=0, color='black', alpha=0.7, linewidth=tick_line_width)

# Remove ticks on x and y axes
plt.xticks([])
plt.yticks([])

plt.show()
