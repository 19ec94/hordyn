import pandas as pd
import glob

files = glob.glob("output/phi_t*.csv")
df = pd.read_csv(files[0])
print("Columns:", list(df.columns))
print("First few rows:")
print(df.head())
print("\nShape:", df.shape)

