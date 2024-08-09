import matplotlib.pyplot as plt

# Define the variables and their properties
variables = {
    'Total Sediment': ('magenta'),
    'River Freshwater Flux': ('green'),
    'Ice Volume': ('red'),
    'Area Cell': ('brown'),
    'River Freshwater Flux NO AC': ('cyan'),
    'Sediment Concentration': ('purple'),
}

# Create a dummy plot to generate the legend
plt.figure(figsize=(10, 2))

# Plot each variable with its corresponding color
for var_name, color in variables.items():
    plt.plot([], [], color=color, label=var_name, linewidth=2)

# Create the legend
plt.legend(loc='center', fontsize=14)

# Hide the axes
plt.gca().set_axis_off()

# Adjust layout
plt.tight_layout()

# Save the legend as an image
plt.savefig('legend_for_powerpoint.png')

# Display the legend
plt.show()

