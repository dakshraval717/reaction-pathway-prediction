import requests
import pandas

cactusDataframe = pandas.read_csv('C:\\Users\\daksh\\Downloads\\Research\\Thermodynamics\\Web Scraping\\Data\\cactus.csv')

cactusSession = requests.Session()

urls = cactusDataframe['Reactant1'].tolist()
reactants1 = []
rowIndex = 0
for url in urls:
    rowIndex += 1
    try:
        response = requests.get(url)
        response.raise_for_status()
        smiles = response.text.strip()
        reactants1.append(smiles)
        print('Successfully converted row', rowIndex, 'to', smiles)
    except:
        print('unsuccessful at row', rowIndex)
        continue