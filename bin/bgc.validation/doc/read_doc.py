import numpy as np
import csv
import datetime
import pickle

# Initialize lists to store the data
dates = []
latitudes = []
longitudes = []
pressures = []
docs = []

# Read and process the CSV file using semicolon as delimiter
with open("doc.csv", newline='', encoding='utf-8') as csvfile:
    reader = csv.reader(csvfile, delimiter=';')
    headers = next(reader)  # Skip header line
    units = next(reader)    # Skip units line

    for row in reader:
        if len(row) < 5:
            continue  # Skip incomplete rows
        try:
            # Replace comma with dot and convert to appropriate types
            date = int(row[0])
            lat = float(row[1].replace(',', '.'))
            lon = float(row[2].replace(',', '.'))
            pressure = float(row[3].replace(',', '.'))
            doc = float(row[4].replace(',', '.'))

            dates.append(date)
            latitudes.append(lat)
            longitudes.append(lon)
            pressures.append(pressure)
            docs.append(doc)
        except ValueError:
            continue  # Skip rows with non-numeric values like 'NA'

# Convert lists to numpy arrays
dates_array = np.array(dates)
latitudes_array = np.array(latitudes)
longitudes_array = np.array(longitudes)
pressures_array = np.array(pressures)
docs_array = np.array(docs)

# Optional: print shapes to verify
print("Dates array shape:", dates_array.shape)
print("Latitudes array shape:", latitudes_array.shape)
print("Longitudes array shape:", longitudes_array.shape)
print("CTD Pressure array shape:", pressures_array.shape)
print("DOC array shape:", docs_array.shape)

index = ~np.logical_or(docs_array < 0.0 ,docs_array > 900.)
doc   = docs_array[index]
lat   = latitudes_array[index] 
lon   = longitudes_array[index]
dates = dates_array[index]
pres  = pressures_array[index]
date=[]
for t in range(dates.shape[0]):
    date.append(datetime.datetime(int(str(dates[t])[0:4]),int(str(dates[t])[4:6]),int(str(dates[t])[6:8])))

f = open('doc.pckl','wb')
pickle.dump([date,lat,lon,pres,doc],f)
f.close()