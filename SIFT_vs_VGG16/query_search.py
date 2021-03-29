import numpy as np
import cv2
from matplotlib import pyplot as plt
import os
import sys
import glob
from timeit import default_timer as timer
from operator import itemgetter

'''
    Usage:
    ./query_search.py -d "database_name" -q "query_imagename" -t "index_type" -r "relevant_images_number"
    
    Example:
    python query_search.py -d COREL -q corel_0000000303_512 -t LINEAR
'''


######## Program parameters
import argparse
parser = argparse.ArgumentParser()

## Database name
parser.add_argument("-d", "--database", dest="db_name",
                    help="input image database", metavar="STRING")

## Query image name
parser.add_argument("-q", "--query", dest="query_name",
                    help="query image name", metavar="STRING")

## Database Index Type
parser.add_argument("-t", "--indextype", dest="indextype",
                    help="index type", metavar="STRING")

## Number of relevant images in the database, considering the query
parser.add_argument("-r", "--relevant", dest="relevant",  type=int,
                    help="relevant image number", metavar="INTEGER", default=4)

args = parser.parse_args()


## Set paths [TO UPDATE]
img_dir="../Images/" + args.db_name + "/"
if(args.db_name == "COREL" or args.db_name == "NISTER"):
	img_dir = "../Images/" + args.db_name + "_queries/"
output_dir="./results/"+ args.db_name
resfilename = "./results/" + args.query_name + "-" + args.indextype


## Load query image
query_filename= img_dir + args.query_name + ".jpg"
if not os.path.isfile(query_filename):
    print "Path to the query "+query_filename+" is not found -- EXIT\n"
    sys.exit(1)

# A utiliser si l'on souhaite utiliser une seule image non pas un dossier d'image en entrée
#queryImage = cv2.imread(query_filename)

#plt.figure(0), plt.title("Image requete")
#plt.imshow(cv2.cvtColor(queryImage, cv2.COLOR_BGR2RGB))

databaseDescriptors = np.load(output_dir + "_DB_Descriptors.npy");
# Chargement de la Map associant numéro d'image et nom d'image de notre base de données
databaseIndex = np.load(output_dir + "_DB_Index.npy");

indexDescriptors = dict(algorithm = 254, filename = output_dir + "_flan_index-LINEAR.dat");
fl=cv2.flann_Index(np.asarray(databaseDescriptors,np.float32),indexDescriptors)

sift = cv2.xfeatures2d.SIFT_create()
# List contenant tous les noms d'image
queryNameList = glob.glob(img_dir+"*.jpg")

mAP = 0
nbquery = min(4,len(queryNameList))

# Parcours de nbquery images d'entrée 
for queryName in queryNameList[:nbquery]:
	
	queryImage = cv2.imread(queryName)
	_, queryDescriptor = sift.detectAndCompute(queryImage,None)

	######## Query search
	start = timer()
	knn=20
	#search_params = dict(checks=50)
	idx, dist=fl.knnSearch(queryDescriptor,knn,params={})

	print "indices \n", idx
	print "distances \n", dist

	vote = []
	# Initialisation d'une map associant nombre d'apparition et path de l'image du descripteur
	for descriptor in databaseIndex:
		if int(descriptor[0]) != len(vote):
			vote.append([0,descriptor[1]])

	# Incrementation des différentes valeurs d'apparition des descripteurs
	for i in idx:
		for j in i:
			vote[int(databaseIndex[j][0])-1][0] = vote[int(databaseIndex[j][0])-1][0] +1
		
	# On trie par ordre décroissant les descripteurs les plus apparus au moins apparus
	sorted_vote = sorted(vote, key=itemgetter(0), reverse = True)

	end = timer()
	print "Search time: " + str(end - start)
	
	########################
	### Display the top images
	########################
	
	top=10
	plt.figure(1), plt.title(args.indextype )
	for i in range(top):
		img = cv2.imread(sorted_vote[i][1])
		score = sorted_vote[i][0]
		plt.subplot(2,5,i+1),plt.imshow(cv2.cvtColor(img, cv2.COLOR_BGR2RGB))
		plt.title('rank '+str(i+1)), plt.xticks([]), plt.yticks([]),plt.xlabel(str(score))

	plt.savefig(resfilename + "_top" + str(top) +".png")
	plt.show()
	
	########################
	### Evaluation
	########################
	
	def getImageId(imname):
		if (args.db_name == "COREL"):
			Id = imname.split('_')[1]
		elif (args.db_name == "NISTER"):
			Id = imname.split('-')[1]
		elif (args.db_name == "Copydays"):
			Id = imname.split('_')[-2]
		else:
			Id = imname.split('.')[-1]
		
		return Id


	base=os.path.basename(queryName)
	queryId=getImageId(os.path.splitext(base)[0])
	print "Identifiant de la requete : ", queryId

	rpFile = open(resfilename + "_rp.dat", 'w')
	precision = np.zeros(top, dtype=float)
	recall = np.zeros(top, dtype=float)

	nbRelevantImage = args.relevant
	nbPertinant = 0
	rang = 1
	API = 0

	# Calcul mAP, précision et rappel
	for i in range(top):
		print recall;
		if(getImageId(sorted_vote[i][1]) == queryId):
			nbPertinant = nbPertinant +1
			rel = 1
		else:
			rel = 0
			
		precision[rang-1] = nbPertinant/float(rang)
		recall[rang-1] = nbPertinant/float(nbRelevantImage)
		rpFile.write(str(precision[rang-1]) + '\t' + str(recall[rang-1]) +  '\n')
		API = API + 1/float(nbRelevantImage)*precision[rang-1]*rel
		
		rang = rang + 1
		
	mAP = mAP + 1/float(nbquery)*API

	# Plot Precision-Recall curve
	plt.clf()
	plt.plot(recall, precision, lw=2, color='navy',
			 label='Precision-Recall curve')
	plt.xlabel('Recall')
	plt.ylabel('Precision')
	plt.ylim([0.0, 1.05])
	plt.xlim([0.0, 1.05])
	plt.title('Precision-Recall for '+args.query_name)
	plt.legend(loc="upper right")
	plt.savefig(resfilename + "_rp.png")
	#plt.savefig(output_dir + args.query_name + "_rp.pdf", format='pdf')
	plt.show()

print mAP


