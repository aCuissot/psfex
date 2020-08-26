import cv2 as cv
import os
from math import log10

"""
On suppose ne pas avoir affaire à des images d'aire > 400000000px²
"""
def int2str(a):
    return str(0)*(3-int(log10(a)))+str(a)

folderIn = './input_dataset'
folderOut = './dataset/'
outputSize = 200


for imageName in os.listdir(folderIn):
    image = cv.imread(folderIN + "/" + imageName)
    cnt = 0
    x_max = image.shape[0]
    y_max = image.shape[1]
    if x_max>=outputSize and y_max>=outputSize:
        for i in range(int(x_max/outputSize)-1):
            for j in range(int(y_max/outputSize)-1):
                cv.imwrite(folderOut + int2str(cnt) + imageName, img[outputSize*i:outputSize*(i+1)][outputSize*j:outputSize*(j+1)])
                cnt += 1
            