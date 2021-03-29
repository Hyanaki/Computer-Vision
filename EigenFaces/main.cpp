#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <visp3/core/vpImage.h>
#include <visp3/io/vpImageIo.h>

#include "EigenFacesDB.cpp"




std::vector<std::string> buildPathImagesAttFaces()
{
	std::vector<std::string> v;
	for(int nbDir=1; nbDir<=40; nbDir++)
		for(int nbImage=1; nbImage<11;nbImage++)
		{
			std::ostringstream ss;
			ss << "../../../att_faces/s" << nbDir << "/" << nbImage << ".pgm";
			v.push_back(ss.str());
		}
	return v;
}

int main()
{
	std::cout << "[INFO] Construction du path ..." << std::endl;
	std::vector<std::string> paths = buildPathImagesAttFaces();
	std::cout << "[INFO] Creation de base de donnees ..." << std::endl;
	EigenFacesDB eigen;

    //eigen.buildBDFaces(paths,400);
	//eigen.writeSynthesisError("../../../att_faces/s11/5.pgm", "../results/Jp.png");




	//Pour les visages tests
	//eigen.buildBDFacesIdentification(paths,400);
	//eigen.writeSynthesisError(paths[eigen.getIdTest()[65]], "../results/Jp.png");



	eigen.distanceOverview(paths,"../results/distanceMatrix.png", 45);
	//eigen.curveFaces(paths,70,3200);
	
	return 0;
}
