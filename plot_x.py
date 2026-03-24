import pandas as pd
import matplotlib.pyplot as plt

# Skip only the duplicate header (1 row after first header)
df = pd.read_csv("pcos/case_log.csv", skiprows=[4])  # Skip the 2nd header row

df_blue = pd.read_csv("Cow419-2ndWave-BlueCurve.csv")  # Skip the 2nd header row
df_orange = pd.read_csv("Cow419-2ndWave-OrangeCurve.csv")  # Skip the 2nd header row
df_red= pd.read_csv("Cow419-2ndWave-RedCurve.csv")  # Skip the 2nd header row
df_green = pd.read_csv("Cow419-2ndWave-GreenCurve.csv")  # Skip the 2nd header row

plt.figure(figsize=(10, 6))
plt.plot(df['t'], df['x0'], color="blue" ,label='x0', linewidth=1.5)
plt.plot(df['t'], df['x1'], color="orange" ,label='x1', linewidth=1.5)
plt.plot(df['t'], df['x2'], color="green" ,label='x2', linewidth=1.5)
plt.plot(df['t'], df['x3'], color="red" ,label='x3', linewidth=1.5)

plt.plot(df_blue['Days']-10, df_blue['Diameter'], color="blue", marker='o', linestyle="", label='x0-exp')
plt.plot(df_orange['Days']-10, df_orange['Diameter'], color="orange", marker='o', linestyle="", label='x1-exp')
plt.plot(df_green['Days']-10, df_green['Diameter'], color="green", marker='o', linestyle="", label='x2-exp')
plt.plot(df_red['Days']-10, df_red['Diameter'], color="red",  marker='o', linestyle="", label='x3-exp')

#plt.plot(df['t'], df['s0'], color="blue" ,label='s0', linestyle='--', linewidth=1.5)
#plt.plot(df['t'], df['s1'], color="orange" ,label='s1', linestyle='--', linewidth=1.5)
#plt.plot(df['t'], df['s2'], color="yellow" ,label='s2', linestyle='--', linewidth=1.5)
#plt.plot(df['t'], df['s3'], color="red" ,label='s3', linestyle='--', linewidth=1.5)

plt.xlabel('t')
plt.ylabel('x')
plt.title('x0, x1, x2, x3 vs t')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

