import pandas
import requests
import regex
import os

cactusDataframe = pandas.read_csv('C:\\Users\\daksh\\Downloads\\Research\\Thermodynamics\\Web Scraping\\Data\\cactus.csv')
fullyProcessedDataframe = cactusDataframe.copy()
semiProcessedDataframe = cactusDataframe.copy()

columns2process = ['Reactant1', 'Reactant2', 'Reactant3', 'Product1', 'Product2', 'Product3']
columns2keep = []

cactusSession = requests.Session()

for rowIndex, rowContents in cactusDataframe.iterrows():

    for columnName in columns2process:
        url = str(rowContents[columnName])
        if url.strip() == '' or not url.strip().startswith('https://cactus.nci.nih.gov/chemical/structure/'):
            continue
        try:
            response = cactusSession.get(url, timeout=15000)
            response.raise_for_status()
            smiles = response.text.strip()

            cactusDataframe.at[rowIndex, columnName] = smiles
            continue
        except:
            fullyProcessedRowBool = False
            continue

    rowIndex += 1
    print(f"Successfully processed row:", {rowIndex})

    if rowIndex % 10 == 0:
        fullyProcessedDataframe = cactusDataframe.copy()
    
    if rowIndex % 20 == 0:
        fullyProcessedDataframe.to_csv('Fully Processed Cactus.csv')
        print('Saved')
        

    # finalDataframe = pd.concat([dataframe, extractedDataframe], axis=1) # this is wrong, extracted has duplicate junk

fullyProcessedDataframe.to_csv('Fully Processed Cactus.csv')
semiProcessedDataframe.to_csv('Semi Processed Cactus')