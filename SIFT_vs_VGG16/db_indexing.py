import os
import sys
import glob
import cv2
import numpy as np
import argparse

from timeit import default_timer as timer

''' 
    Usage :
    ./db_indexing.py -d "database_name"
    
    Example :
    ./db_indexing.py -d "base1"
'''

######## Program parameters

parser = argparse.ArgumentParser()

## Database name
parser.add_argument("-d", "--database", dest="db_name",
                    help="input image database", metavar="STRING", default="None")

args = parser.parse_args()

## Set paths
img_dir="./../Images/" + args.db_name + "/"
imagesNameList = glob.glob(img_dir+"*.jpg")
output_dir="./results/" + args.db_name

if not os.path.exists(img_dir):
    print "The directory containing images: "+img_dir+" is not found -- EXIT\n"
    sys.exit(1)


start = timer()
databaseDescriptors = []
databaseIndex = []
nbimages = 0
nbdescriptors = 0
for imageName in imagesNameList:
	img = cv2.imread(imageName)
	sift = cv2.xfeatures2d.SIFT_create()
	_, des = sift.detectAndCompute(img,None)			
	nbimages = nbimages + 1
	if (des is not None):
		for descriptor in des:
			databaseDescriptors.append(descriptor)
			databaseIndex.append([nbimages,imageName])
			nbdescriptors = nbdescriptors + 1
			
print databaseIndex[1]
FLANN_INDEX_ALGO=0
index_params = dict(algorithm = FLANN_INDEX_ALGO)   # for linear search
### OpenCV 2.4.13
print index_params
fl=cv2.flann_Index(np.asarray(databaseDescriptors,np.float32),index_params)

end = timer()
print "Indexing time: " + str(end - start)

np.save(output_dir + "_DB_Descriptors.npy", databaseDescriptors)
np.save(output_dir + "_DB_Index.npy", databaseIndex)
print "Nb Images : " + str(nbimages) + ", Nb Descriptors : "+ str(nbdescriptors)

fl.save(output_dir + "_flan_index-LINEAR.dat")
 
