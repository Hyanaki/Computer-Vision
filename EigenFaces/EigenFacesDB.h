/*
 * EigenFacesDB.h
 *
 *  Created on: Dec 17, 2010
 *      Author: beltegeuse
 */

#ifndef EIGENFACESDB_H_
#define EIGENFACESDB_H_

#include <visp/vpColVector.h>
#include <visp/vpMatrix.h>
#include <vector>
#include <string>

class EigenFacesDB
{
private:
	/**
	 * Attributs
	 */
	// Vecteur qui represente notre image moyenne
	vpColVector m_vMean;
	// Matrice contenant les vecteur U
	vpMatrix m_U;
	// Vecteur contenant les valeurs propres
	vpColVector m_eigenValue;
	// Les caracteristiques de l'image
	int m_h, m_w;
	int m_size;
	// Le nombre de vecteur U
	int m_maxEigenFace;

	vpColVector distanceMatrix;
	std::vector<vpRowVector> imagesDbRef;
    std::vector<vpRowVector> imagesDbTest;

	std::vector<vpColVector> imagesMeanMatrices;

	std::vector<int> idRef;
    std::vector<int> idTest;

public:
	/**
	 *  Constructeurs et destructeurs
	 */
	EigenFacesDB();
	virtual ~EigenFacesDB();

	// Construire notre base de donnee
	void buildBDFaces(const std::vector<std::string>& paths, int maxEigenFace = 50);
	// Ecrire des images sur le systeme
	// * L'image moyenne
	void writeMeanImage(const std::string& path) const;
	// * les images des eigenFaces
	void writeEigenFacesImage(const std::string& directory = "");
	// * l'image synthetisee
	vpColVector computeCoordonateImage(const vpMatrix &U, const vpColVector &A, int nbDim);
	void writeSynthesisImage(const std::string& image, const std::string& out, int nbDim, vpColVector &Jp, vpColVector &J);
	void writeSynthesisError(const std::string& image, const std::string& path);
	double computeSynthesisError(const std::string& image, int K);
	// Calculer le vecteur pour identifier une image
	vpColVector computeIdVector(const std::string& image, int nbDim);
	// Pour ecrire les valeurs propres
	void writeValeursPropres(const vpMatrix &A);
	// Matrice de distance
	double computeDistanceImage(int k, int i, int j);
	double distanceIdVectors(vpColVector& id1, vpColVector& id2);
	void writeMatriceDistance(std::vector<std::string> paths, const std::string& out, int nbDim, int k);
	//compute min and max distance intra and inter classes
	void distanceOverview(std::vector<std::string> paths, const std::string& out, int nbDim);
	void curveFaces(std::vector<std::string> paths, int nbDim,int seuil);
	void buildBDFacesIdentification(const std::vector<std::string>& paths, int maxEigenFace);

	std::vector<int> getIdTest(){
		return idTest;
	}
private:
	void buildMeanImage(const std::vector<vpRowVector>& imagesDb, vpColVector& averageFace);
	vpMatrix buildMatrixA(const std::vector<vpRowVector> &imagesDb);
	void computeU(vpMatrix& A, vpMatrix& eigenVec);
	void writeEigenFaceInImage(int m_w, int m_h,const std::string& path, vpColVector &v);
	void writeFaceCenter(vpMatrix &A, int nbFaces, int w, int h);
	void computeCoordonateImage(const vpMatrix &U, const vpMatrix &A );
	
};
#endif /* EIGENFACESDB_H_ */
