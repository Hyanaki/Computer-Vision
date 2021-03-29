#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <visp3/core/vpImage.h>
#include <visp3/io/vpImageIo.h>
#include "EigenFacesDB.h"


#define VP_TRACE

//! \example tutorial-viewer.cpp
//! [Include display]
#include <visp3/gui/vpDisplayD3D.h>
#include <visp3/gui/vpDisplayGDI.h>
#include <visp3/gui/vpDisplayGTK.h>
#include <visp3/gui/vpDisplayX.h>
#include <visp3/gui/vpDisplayOpenCV.h>
#include <visp3/core/vpPoint.h>
#include <visp3/core/vpPlane.h>
#include <visp3/gui/vpPlot.h>
#include <visp3/core/vpMeterPixelConversion.h>

#include <limits>
#include <stdlib.h>

EigenFacesDB::EigenFacesDB()
{

}

EigenFacesDB::~EigenFacesDB()
{
    
}

void  EigenFacesDB::writeEigenFaceInImage(int m_w, int m_h, const std::string& path, vpColVector &v)
{
    double maxV = v.getMaxValue();
    double minV = v.getMinValue();
    vpImage<unsigned char> I(m_h,m_w);
    for(int i=0; i<m_h; i++){
        for(int j = 0; j<m_w; j++)
        {
            I[i][j] = ((v[i*m_w+j]-minV)/(maxV-minV))*255;
        }
    }
    vpImageIo::write(I, path);
}

void EigenFacesDB::buildBDFaces(const std::vector<std::string>& paths, int maxEigenFace)
{
    int m_maxEigenFace = std::min(maxEigenFace, (int)paths.size());
    
    // On calcul les attibuts de l'image
    vpImage<unsigned char> I;
    vpImageIo::read(I,*paths.begin());
    m_w = I.getWidth();
    m_h = I.getHeight();
    m_size = m_h*m_w;
    
    std::cout << " * Caracteristique de l'images : " << m_h << "x" << m_w << std::endl;
    std::cout << " * Nombre d'image de la base : " << paths.size()<< std::endl;
    std::cout << " * Nombre de U : " << m_maxEigenFace << std::endl;
    
    std::vector<vpRowVector> imagesDb(m_maxEigenFace);
    for(int indice = 0;indice < m_maxEigenFace ; indice++)
    {
        vpImageIo::read(I,paths[indice]);
        imagesDb[indice] = vpRowVector(m_w*m_h);
        for(int j =0;j<m_w;j++){
            for(int i =0;i<m_h;i++){
               imagesDb[indice][j + m_w*i] = (double)I[i][j];
            }
        }
    }
    
    // Creation du vpColVector pour le mean face
    std::cout << "[INFO] Creation Mean images ...." << std::endl;
    m_vMean = vpColVector(m_size);
    buildMeanImage(imagesDb,m_vMean);
    writeEigenFaceInImage(m_w,m_h,"../results/meanFace.png", m_vMean);
    // Calcul de la matrice A
    std::cout << "[INFO] Calcul de A ... " << std::endl;
    vpMatrix A = buildMatrixA(imagesDb);
    
    writeFaceCenter(A,10, m_w, m_h);
    std::cout << "[INFO] Fin calcul BD ... " << std::endl;
    vpColVector v(m_w*m_h);
    writeValeursPropres(A);
    std::string path;
    for (int i=0; i<10;i++){
        std::stringstream ss;
        ss << "../results/EigenFaces" << i << ".png";
        path = ss.str();
        v = m_U.getCol(i);
        writeEigenFaceInImage(m_w,m_h, path, v);
    }

    // Créer le VpColVector A (J-Mean) pour toutes les images du path 
    for(int indice = 0;indice < m_maxEigenFace ; indice++)
    {
        vpImageIo::read(I,paths[indice]);
        vpColVector A1(m_w*m_h);
        for(int j =0;j<m_w;j++){
            for(int i =0;i<m_h;i++){
                A1[j + m_w*i] = (double)I[i][j] - m_vMean[j + m_w*i];
            }
        }
        imagesMeanMatrices.push_back(A1);
    }
}

void EigenFacesDB::buildMeanImage(const std::vector<vpRowVector>& imagesDb, vpColVector& averageFace){
    for (int i = 0; i<imagesDb.size(); i++){
        for (int j=0; j<imagesDb[0].size(); j++){
            averageFace[j] = averageFace[j]+(double)imagesDb[i][j]/imagesDb.size();
        }
    }
}

vpMatrix EigenFacesDB::buildMatrixA(const std::vector<vpRowVector> &imagesDb){
    
    vpMatrix A(imagesDb[0].size(),imagesDb.size());
    for (int i = 0; i<imagesDb.size(); i++){
        for (int j=0; j<imagesDb[0].size(); j++){
            A[j][i] = (double)imagesDb[i][j]-m_vMean[j];
        }
    }
    return A;
}

void EigenFacesDB::writeFaceCenter(vpMatrix &A, int nbFaces, int w, int h){
    
    vpColVector result(A.getRows());
    //vpColVector difference(A.getCols());
    std::string path = "";
    for (int i=0; i<nbFaces; i++){
        std::stringstream ss;
        for (int j = 0; j<A.getRows(); j++){
            result[j] = A[j][i];
        }
        ss << "../results/A" << i << ".png";
        path = ss.str();
        std::cout << path << std::endl;
        writeEigenFaceInImage(w,h,path,result);
        ss.clear();
    }

}

void EigenFacesDB::writeValeursPropres(const vpMatrix &A){   
    // On calcul les attibuts de l'image
    vpColVector W;
    m_U= A;
    vpMatrix V;
    m_U.svd(W,V);
    m_eigenValue= vpColVector(W.size());
    for (int i = 0; i<W.size();i++){
        m_eigenValue[i]=W[i]*W[i];
        if (i<10){
            std::cout<<m_eigenValue[i]<<std::endl;;
        }
    }
    
}

vpColVector EigenFacesDB::computeCoordonateImage(const vpMatrix &U, const vpColVector &A, int nbDim){
    vpColVector W(nbDim);
    for (int i=0; i<nbDim; i++){
        for (int j = 0; j<U.getRows(); j++){
            W[i] +=  U[j][i]*A[j];
        }
    } 
    return W;   
     
}

void EigenFacesDB::writeSynthesisImage(const std::string& image, const std::string& out, int nbDim, vpColVector &Jp, vpColVector &Jmean){
    Jp = m_vMean;
    vpColVector W= computeCoordonateImage(m_U,Jmean,nbDim);
    for (int k=0; k<nbDim; k++){
        Jp += W[k]*m_U.getCol(k);
    }
    writeEigenFaceInImage(m_w,m_h,out, Jp);
    
}
void EigenFacesDB::writeSynthesisError(const std::string& image, const std::string& path){
    vpPlot plot(4, 700, 700, 100, 200, "Curves...");
    vpImage<unsigned char> I;
    vpImageIo::read(I,image);
    vpDisplayX d;
    d.init(I);
    vpDisplay::display(I);
    vpColVector Jmean(m_w*m_h);
    vpColVector Jp(m_w*m_h);

    for(int i =0;i<m_h;i++){
        for(int j =0;j<m_w;j++){
           Jmean[j + m_w*i] = (double)I[i][j] - m_vMean[j + m_w*i];
        }
    }

    char title[40];
    strncpy( title, "Eigen values", 40 );
    plot.setTitle(0,title);
    plot.initGraph(0,1);

    strncpy( title, "Reconstruction Error", 40 );
    plot.setTitle(1, title);
    plot.initGraph(1,1);

    double cumulSum = 0;
    m_eigenValue = m_eigenValue/(double)m_eigenValue.sum();
    
    vpColVector e(Jmean.size());
    for (int i = 0; i<m_eigenValue.size(); i++){
        cumulSum += m_eigenValue[i];
        plot.plot(0,0,i, cumulSum) ;
        writeSynthesisImage(image, path,i, Jp,Jmean);
        e = Jmean+m_vMean-Jp;
        plot.plot(1,0,i, e.euclideanNorm()) ;
    }
    vpColVector J = Jmean+m_vMean;
    vpDisplay::getClick(I);
    writeEigenFaceInImage(m_w,m_h,"../results/J.png", J);
    plot.saveData(0,"cumulValue.txt","#");
    plot.saveData(1,"Error.txt","#");
    std::cout << "Clicker sur l'image pour terminer" << std::endl ;

    
    
}

double EigenFacesDB::computeDistanceImage(int k, int i, int j){

    vpColVector W1 = computeCoordonateImage(m_U,imagesMeanMatrices[i],k);
    vpColVector W2 = computeCoordonateImage(m_U,imagesMeanMatrices[j],k);
    return (W2-W1).euclideanNorm();

};
//double distanceIdVectors(vpColVector& id1, vpColVector& id2);
void EigenFacesDB::writeMatriceDistance(std::vector<std::string> paths, const std::string& out, int nbDim, int K_eigen){
    
    double distance;
    distanceMatrix= vpColVector(nbDim*nbDim);
    
    for (int i =0; i <nbDim; i++){
        for (int k = i; k<nbDim; k++){
            distance = computeDistanceImage(K_eigen,i,k);
            distanceMatrix[i*nbDim+k] = distance; 
            distanceMatrix[k*nbDim+i] = distance;          
        }
    }
    
    writeEigenFaceInImage(nbDim, nbDim, out, distanceMatrix);
    

}

void EigenFacesDB::distanceOverview(std::vector<std::string> paths, const std::string& out, int nbDim){
 //To answer Question 13-14

    buildBDFaces(paths,nbDim);
    double distanceInter, distanceIntra;
    double dminInter= std::numeric_limits<double>::max();
    double dmaxInter=std::numeric_limits<double>::min();
    double dminSame=std::numeric_limits<double>::max();
    double dmaxSame=std::numeric_limits<double>::min();
    writeMatriceDistance(paths, out,nbDim,20);

    for (int i =0; i <nbDim; i++){
        for(int k=i+1; k<nbDim; k++){
            if (k<i+10){
                distanceIntra = distanceMatrix[i*nbDim +k];
                if( distanceIntra<dminSame){
                    dminSame = distanceIntra;
                }
                else if(distanceIntra>dmaxSame){
                    dmaxSame = distanceIntra;
                }
            }
            else {
                distanceInter = distanceMatrix[i*nbDim+k];
                if (distanceInter<dminInter){
                    dminInter = distanceInter;
                }   
                else if(distanceInter>dmaxInter){
                    dmaxInter = distanceInter;
                } 
            }

        }
        
    }


//To answer Question 13-14
   std::cout<<"The minimum inter class distance is equal to "<<dminInter<<" and the maximum to "<<dmaxInter<<std::endl;
   std::cout<<"The minimum intra class distance is equal to "<<dminSame<<" and the maximum to "<<dmaxSame<<std::endl;
}

void EigenFacesDB::curveFaces(std::vector<std::string> paths, int nbDim,int seuil){
    
    buildBDFacesIdentification(paths,nbDim);
    vpImage<unsigned char> I;
    
    vpImageIo::read(I,paths[0]);

    vpDisplayX d;
    d.init(I);
    vpDisplay::display(I);

    int known;
    int testSize = imagesDbTest.size(); 
    int refSize = imagesDbRef.size();
    double distance;

    vpPlot plot(1, 700, 700, 100, 200, "Curves...");
    char title[40];
    strncpy( title, "Pourcentage de visage reconnu", 40 );
    plot.setTitle(0,title);
    plot.initGraph(0,1);
   // std::cout<< m_U.size()<<std::endl;
    for(int k = 1; k< m_U.size()/m_size;k++){
        known = 0;
        for (int i =0; i <testSize; i++){
            for (int j = 0; j<refSize; j++){
                distance = computeDistanceImage(k,idTest[i],idRef[j]);
                if(distance<seuil){
                    known+=1;
                }    
            }
        }
        plot.plot(0,0,k,(double)known/testSize) ;
    }
    vpDisplay::getClick(I);
}

void EigenFacesDB::buildBDFacesIdentification(const std::vector<std::string>& paths, int maxEigenFace)
{
    int m_maxEigenFace = std::min(maxEigenFace, (int)paths.size());
    srand(time(NULL));
    // On calcul les attibuts de l'image
    vpImage<unsigned char> I;
    vpImageIo::read(I,*paths.begin());
    m_w = I.getWidth();
    m_h = I.getHeight();
    m_size = m_h*m_w;
    int r;
    

    std::cout << " * Caracteristique de l'images : " << m_h << "x" << m_w << std::endl;
    std::cout << " * Nombre d'image de la base : " << m_maxEigenFace<< std::endl;
    

    vpRowVector temp =vpRowVector(m_w*m_h);
    for(int indice = 0;indice < m_maxEigenFace ; indice++)
    {
        vpImageIo::read(I,paths[indice]);
        r = rand() % 100;
        if(r<35){
            for(int j =0;j<m_w;j++){
                for(int i =0;i<m_h;i++){
                    temp[j + m_w*i] = (double)I[i][j];
                }
            }
            imagesDbTest.push_back(temp);
            idTest.push_back(indice);
        }
        else{
            for(int j =0;j<m_w;j++){
                for(int i =0;i<m_h;i++){
                    temp[j + m_w*i] = (double)I[i][j];
                }
            }
             imagesDbRef.push_back(temp);
             idRef.push_back(indice);
             
        }
        
    }
    
    std::cout<<"Database Size : "<< imagesDbRef.size() <<std::endl;
    std::cout<<"Database Size for the Test : "<< imagesDbTest.size() <<std::endl;
    
    // Creation du vpColVector pour le mean face
    std::cout << "[INFO] Creation Mean images ...." << std::endl;
    m_vMean = vpColVector(m_size);
    buildMeanImage(imagesDbRef,m_vMean);
    writeEigenFaceInImage(m_w,m_h,"../results/meanFace.png", m_vMean);
    // Calcul de la matrice A
    std::cout << "[INFO] Calcul de A ... " << std::endl;
    vpMatrix A = buildMatrixA(imagesDbRef);
    
    std::cout << "[INFO] Fin calcul BD ... " << std::endl;

    writeValeursPropres(A);
    std::cout<<"U size checkup : "<< m_U.size()/m_size <<std::endl;


    // Créer le VpColVector A (J-Mean) pour toutes les images du path 
    for(int indice = 0;indice < m_maxEigenFace ; indice++)
    {
        vpImageIo::read(I,paths[indice]);
        vpColVector A1(m_w*m_h);
        for(int j =0;j<m_w;j++){
            for(int i =0;i<m_h;i++){
                A1[j + m_w*i] = (double)I[i][j] - m_vMean[j + m_w*i];
            }
        }
        imagesMeanMatrices.push_back(A1);
    }
}

    