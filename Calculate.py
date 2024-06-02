
###########################calculate LOSH#######################################
import clusterpy
import numpy as np
import shapefile
import pandas as pd

# Function to calculate Local Spatial Heteroscedasticity (LOSH)
def calculateLOSH(keyList, keyDict, a, dataDictionary, dataLength):
    # Calculate E1 values
    E1 = []
    for x in range(dataLength):
        keyList1 = keyDict[x]
        sum = 0
        for i in keyList1:
            sum += dataDictionary[i]
        neighborNumber = len(keyList1)
        localMean = float(sum) / neighborNumber
        selfresidual = dataDictionary[x] - localMean
        E1.append(selfresidual)

    # Calculate H_one value
    h1 = 0
    for x in E1:
        h1 += x ** 2
    H_one = float(h1) / dataLength

    # Calculate LOSH value for a given keyList
    count = 0
    for i in keyList:
        localresidual = E1[i]
        Esq = localresidual ** a
        count += Esq
    losh1 = count
    neighborNumber = len(keyList)
    LOSH = round(float(losh1) / (H_one * neighborNumber), 9)

    return LOSH

# Load data
# Import shapefile data for the regions
Rwanda = clusterpy.importArcData(".../1KM_2013/Rwan_poly_2013")
#Rwanda = clusterpy.importArcData(".../2014/Kigali_polygon_1km_2014")
wt = Rwanda.Wqueen

# Load shapefile for region boundaries
shape = shapefile.Reader(".../Rwanda/1KM_2013/Rwan_poly_2013")
#shape = shapefile.Reader(".../2014/Kigali_polygon_1km_2014")

# Load population data
data = pd.read_csv(".../1KM_2013/ppp_RWA_2013_1km_Aggregated_UNadj.csv")  ##pop = data.Z
#data = pd.read_csv(".../2014/Kigali_2014_1km.csv") ##pop = data.r__2014
id = data.FID
pop = data.Z

# Create a dictionary mapping region IDs to population values
dataDictionary = {}
for a, b in enumerate(id):
    b = int(b)
    dataDictionary[b] = pop[a]

dataLength = len(shape)

# Calculate LOSH for each region
result = []
for x in range(dataLength):
    Keylist = wt[x]
    losh = calculateLOSH(Keylist, wt, 2, dataDictionary, dataLength)
    result.append(losh)
#print(result)

# Significance test using Monte Carlo permutation
areaKeys = dataDictionary.keys()
plist = []
for x in areaKeys:
    Nlist = list(range(0, dataLength))
    betterClusters = 0
    number = len(wt[x])
    Nlist.pop(x)
    for j in range(99):
        permKey = np.random.choice(Nlist, number, False)
        permKey = permKey.tolist()
        randomH = calculateLOSH(permKey, wt, 2, dataDictionary, dataLength)
        if result[x] < randomH:
            betterClusters += 1
    pValue = (betterClusters + 1) / 100.00
    plist.append(pValue)

# Create a DataFrame to store results and p-values
df_LOSH = pd.DataFrame()
df_LOSH['Hi'] = result
df_LOSH['P_sim'] = plist
#df_LOSH.to_csv(".../LOSH_results_2013.csv")
print("LOSH finished")


###########################calculate Local G#######################################
import clusterpy
import numpy as np
import shapefile
import pandas as pd

# Function to calculate the local Getis-Ord G* statistic
def calculateGetisG(keyList, dataMean, dataStd, dataDictionary, dataLength):
    """
    This function returns the local G statistic for a given region.
    keyList is the list of keys of neighbors
    dataLength is the total number of input data units
    dataMean = global mean
    dataStd = global standard deviation
    """
    sum = 0
    for i in keyList:
        sum += dataDictionary[i]
    neighborNumber = len(keyList)
    numerator = sum - (dataMean * neighborNumber)
    denominator = dataStd * ((float(dataLength * neighborNumber - (neighborNumber ** 2)) / (dataLength - 1)) ** 0.5)
    G = numerator / denominator
    return G

# Load data
# Import shapefile data for the regions
rwanda = clusterpy.importArcData(".../1KM_2013/Rwan_poly_2013")
wt = rwanda.Wqueen

# Create a new weight matrix that includes each region itself
new_wt = {}
list1 = []
list2 = []
for x in wt:
    wtlist = wt[x]
    wtlist.append(x)
    list1.append(wtlist)
    list2.append(x)
for a, b in enumerate(list2):
    b = int(b)
    new_wt[b] = list1[a]

# Get data dictionary
shape = shapefile.Reader(".../1KM_2013/Rwan_poly_2013")
data = pd.read_csv(".../1KM_2013/ppp_RWA_2013_1km_Aggregated_UNadj.csv")
ori_data = pd.read_csv(".../1KM_2004/ppp_RWA_2004_1km_Aggregated_UNadj.csv")
id = data.FID
pop = data.X
oripop = ori_data.X  # Set the 1st year population data
ori_Dic = {}
dataDictionary = {}
for a, b in enumerate(id):
    b = int(b)
    dataDictionary[b] = pop[a]
    ori_Dic[b] = oripop[a]

areaKeys = dataDictionary.keys()

# Calculate data mean and standard deviation
dataMean = np.mean(np.double(list(dataDictionary.values())))
dataStd = np.std(np.double(list(dataDictionary.values())))

# Get data length
dataLength = len(shape)

# Calculate first year data mean and standard deviation
fst_dataMean = np.mean(np.double(list(ori_Dic.values())))
fst_dataStd = np.std(np.double(list(ori_Dic.values())))

# Calculate local G* values for each region
clusterGstrValues = {}
resultstr = []
for x in range(dataLength):
    keyList = new_wt[x]
    currentG = calculateGetisG(keyList, fst_dataMean, fst_dataStd, dataDictionary, dataLength)
    resultstr.append(currentG)
    clusterGstrValues[x] = currentG

# Monte Carlo permutation test for significance
plist = []
for x in areaKeys:
    Nlist = list(range(0, dataLength))
    betterClusters = 0
    number = len(new_wt[x][:-1])
    Nlist.pop(x)
    for j in range(999):
        permKey = np.random.choice(Nlist, number, False)
        permKey = permKey.tolist()
        permKey.append(x)
        randomG = calculateGetisG(permKey, fst_dataMean, fst_dataStd, ori_Dic, dataLength)
        if clusterGstrValues[x] >= 0:
            if clusterGstrValues[x] < randomG:
                betterClusters += 1
        else:
            if clusterGstrValues[x] > randomG:
                betterClusters += 1
    pValue = (betterClusters + 1) / 1000.00
    plist.append(pValue)

# Save results to a DataFrame and CSV file
df_G = pd.DataFrame()
df_G['G_str'] = resultstr
df_G['P_sim'] = plist
#df_G.to_csv(".../G_results_2013.csv")
print("Local G finished")

###########################calculate UDI#######################################
# Create UDI column based on G and LOSH results
UDI = []
for i in range(len(df_G)):
    G_value = df_G['G_str'][i]
    LOSH_value = df_LOSH['Hi'][i]
    G_significance = df_G['P_sim'][i] <= 0.05
    LOSH_significance = df_LOSH['P_sim'][i] <= 0.05
    LOSH_category = 'Not Significant'

    if LOSH_significance:
        if LOSH_value > 1:
            LOSH_category = 'High'
        elif LOSH_value < 1:
            LOSH_category = 'Low'
    else:
        LOSH_category = 'Not Significant'

    if G_value < 0:
        if LOSH_category == 'Low':
            UDI.append(1)
        elif LOSH_category == 'Not Significant':
            UDI.append(2)
        else:
            UDI.append(3)
    elif G_value > 0:
        if LOSH_category == 'High':
            UDI.append(7)
        elif LOSH_category == 'Not Significant':
            UDI.append(8)
        else:
            UDI.append(9)
    else:
        if LOSH_category == 'High':
            UDI.append(4)
        elif LOSH_category == 'Not Significant':
            UDI.append(5)
        else:
            UDI.append(6)

# Create DataFrame to store results
df_result = pd.DataFrame()
df_result['G_str'] = df_G['G_str']
df_result['G_P_sim'] = df_G['P_sim']
df_result['LOSH'] = df_LOSH['Hi']
df_result['LOSH_P_sim'] = df_LOSH['P_sim']
df_result['UDI'] = UDI

# Save DataFrame to CSV
df_result.to_csv(".../UDI_results_2013.csv", index=False)

print("Analysis finished")
