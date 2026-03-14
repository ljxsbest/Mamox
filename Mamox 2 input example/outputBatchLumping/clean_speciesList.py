import pandas as pd

# Read species to remove from text file (split by whitespace)
species_to_remove = set()
with open("species_clean_Mech25_37_core.txt", "r", encoding="utf-8") as f:
    for line in f:
        tokens = line.strip().split()
        species_to_remove.update(tokens)

# Read the CSV file
df = pd.read_csv("speciesList.csv")

# Check if the column exists
if "Mech. name" not in df.columns:
    raise ValueError("Column 'Mech. name' not found in the CSV file. Please check the column name.")

# Filter out rows where "Mech. name" is in the removal list
df_clean = df[~df["Mech. name"].isin(species_to_remove)]

# Save to a new CSV file
df_clean.to_csv("speciesList_clean.csv", index=False, encoding="utf-8")

print("Done. The cleaned file has been saved as speciesList_clean.csv")
