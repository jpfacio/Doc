import pandas as pd
import json

with open('Data/Chen_Analysis/chen_data.json', 'r') as f:
    data = json.load(f)
    
records = data['data']

df = pd.DataFrame(records)
df.to_csv('Data/Chen_Analysis/chen_data.csv', index=False)