import pandas as pd

df = pd.read_csv("Data/Chen_Analysis/chen_data.csv")
df = df.drop(["sample_type", "environment", "host", "technology", "country", "city",
              "collected_date", "isolation_source", "platform", "altitude", "host_disease"], axis=1)
df = df[df["latitude_and_longitude"] != "missing"]
filtered = df.sample(frac=0.05, random_state=5)

print(f"Number of lines in raw data: {df.shape[0]}")
print(f"Number of lines in filtered data: {filtered.shape[0]}")

filtered.to_csv("Data/Chen_Analysis/chen_data_5.csv", index=False)

N_s = filtered["latitude_and_longitude"].str.count("S").sum()
print(f"Samples from South Hemisphere: {N_s}")

filtered["download"] = filtered.apply(
    lambda x: f"{x['download']}{x['assembly_id']}/{x['sample_name']}.fa.gz", 
    axis=1
)

links = filtered["download"].tolist()

wget = [f"wget -P Data/Chen_Analysis/Bins -c {link}" for link in links]

with open("Scripts/Chen_Analysis/links.sh", "w") as f:
    for cmd in wget:
        f.write(cmd + "\n")

