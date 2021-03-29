
from keras.applications.vgg16 import VGG16
from keras.preprocessing import image
from keras.applications.vgg16 import preprocess_input
from keras.models import Model
import numpy as np

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

model = VGG16(weights='imagenet', include_top=False, pooling='max')

start = timer()
databaseDescriptors = []
databaseIndex = []
nbimages = 0
for imageName in imagesNameList:
	img = image.load_img(imageName, target_size=(224, 224))
	x = image.img_to_array(img)
	
	x = np.expand_dims(x, axis=0)
	x = preprocess_input(x) 
	
	des = model.predict(x)
	
	databaseDescriptors.append(des[0])
	databaseIndex.append([nbimages,imageName])
	
	nbimages = nbimages + 1
			
print databaseIndex[1]
FLANN_INDEX_ALGO=0
index_params = dict(algorithm = FLANN_INDEX_ALGO)   # for linear search
### OpenCV 2.4.13
print index_params
fl=cv2.flann_Index(np.asarray(databaseDescriptors,np.float32),index_params)

end = timer()
print "Indexing time: " + str(end - start)

np.save(output_dir + "_DB_VGG16_Descriptors.npy", databaseDescriptors)
np.save(output_dir + "_DB_VGG16_Index.npy", databaseIndex)
print "Nb Images : " + str(nbimages) + ", Nb Descriptors : "+ str(nbimages)

fl.save(output_dir + "_VGG16_flan_index-LINEAR.dat")

