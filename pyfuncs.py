# parse coordinate strings to ndarry(?)
import numpy as np
import pandas as pd
def parse_coordinate(coordinateString, axisLim=5):
    splitted = coordinateString.split(";")
    splitted.pop() # discard the last element because its empty
    result = np.zeros((axisLim, axisLim))
    for i in splitted:
        secondSplit = i.split(",")
        result[int(secondSplit[0]), int(secondSplit[1])] += 1
    return result

# plot out coordinate arrays
import matplotlib.pyplot as plt
import pandas as pd
def coordinate_plot(csvFile="yeast_test_parsing_df.csv", plotNum=5,
        plotSizeColor=5):
    df = pd.read_csv(csvFile)
    temp = df["coordinates"].sample(plotNum)
    print(temp)
    # temp = df["coordinates"].loc[[2731, 2783, 2730, 2784, 2740]]
    temp.index=range(plotNum)
    result = [parse_coordinate(x) for x in temp]
    for i in range(plotNum):
        plt.subplot(plotNum, plotNum, i+1)
        plt.xlim(-1,plotSizeColor)
        plt.ylim(-1,plotSizeColor)
        print(result[i])
        # plt.imshow(result[i], vmin=0, vmax=plotSizeColor)
        plt.imshow(result[i])
    plt.show()

# temporary code
# ==============
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
df = pd.read_csv("./yeast_for_tsne.csv")
temp = df
temp.index = df["name"]
temp = temp.drop("name", axis=1)
tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=5000)
temp1 = tsne.fit_transform(temp)
temp2 = pd.DataFrame(temp1)
temp2.index = df.index
# temp1 is for visualization
# temp2 contains all genes names as index names
plt.close("all")
plt.scatter(temp1[:,0], temp1[:,1])
plt.show(block=False)
temp2.columns = ["X", "Y"]
print(temp2[temp2.index.str.contains("Rad53")])
print(temp2[temp2.index.str.contains("Rec114")])
print(temp2[temp2.index.str.contains("Rad9")])
print(temp2[temp2.index.str.contains("Chk1")])
print(temp2[temp2.index.str.contains("Rec8")])
print(temp2[temp2.index.str.contains("Rad51")])
print(temp2[temp2.index.str.contains("Sin3")])
