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

# temporary code
# plot out coordinate arrays
import matplotlib.pyplot as plt
import pandas as pd
def coordinate_plot(csvFile="yeast_test_parsing_df.csv", plotNum=5,
        plotSizeColor=5):
    df = pd.read_csv(csvFile)
    # temp = df["coordinates"].sample(plotNum)
    temp = df["coordinates"].loc[[2731, 2783, 2730, 2784, 2740]]
    temp.index=range(plotNum)
    result = [parse_coordinate(x) for x in temp]
    for i in range(plotNum):
        plt.subplot(plotNum, plotNum, i+1)
        plt.xlim(-1,plotSizeColor)
        plt.ylim(-1,plotSizeColor)
        print(result[i])
        plt.imshow(result[i], vmin=0, vmax=plotSizeColor)
    plt.show()

