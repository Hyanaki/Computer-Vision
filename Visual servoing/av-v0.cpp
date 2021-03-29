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
#include <visp3/core/vpCameraParameters.h>

#include <visp3/core/vpExponentialMap.h>
#include <visp3/core/vpVelocityTwistMatrix.h>



//! [Include display]
//! [Include io]
#include <visp3/io/vpImageIo.h>
//! [Include io]

using namespace std ;


void
display(vpCameraParameters& cam, vpImage<unsigned char> &I, vpColVector &x, vpColVector &xd)
{
    for (int i = 0 ; i < x.getRows()/2 ; i++)
    {
        vpImagePoint u,ud ;
        vpMeterPixelConversion::convertPoint(cam,x[2*i],x[2*i+1],u) ;
        vpDisplay::displayPoint(I,u,vpColor::red, 2) ;

        vpMeterPixelConversion::convertPoint(cam,xd[2*i],xd[2*i+1],ud) ;
        vpDisplay::displayCircle(I,ud,10,vpColor::green) ;

    }

    vpDisplay::flush(I) ;
}


// Projection d'un point 3D sur le plan image  X(3), x(2)
void
project(vpColVector &X, vpColVector &x)
{
	x[0] = X[0]/X[2];
	x[1] = X[1]/X[2];
}

// Changement de repere bX(3), aX(3), aTb est une matrice homogène
void changeFrame(const vpColVector &bX, const vpHomogeneousMatrix &aTb,  vpColVector &aX)
{
	vpColVector bX4(4);
	bX4[0] = bX[0];
	bX4[1] = bX[1];
	bX4[2] = bX[2];
	bX4[3] = 1;
	aX = aTb*bX4;
}

// Calcul de la matrice d'interaction de point 2D
void
computeInteractionMatrix(vpColVector * cX, vpColVector * xd, vpMatrix &Lx, int nbPoints)
{
    for (int i=0;i<nbPoints;i++){
            Lx[2*i][0] = -1/cX[i][2];
            Lx[2*i][1] = 0;
            Lx[2*i][2] = xd[i][0]/cX[i][2];
            Lx[2*i][3] = xd[i][0]*xd[i][1];
            Lx[2*i][4] = -(1+pow(xd[i][0],2));
            Lx[2*i][5] = xd[i][1];
            
            Lx[2*i+1][0] = 0;
            Lx[2*i+1][1] = -1/cX[i][2];
            Lx[2*i+1][2] = xd[i][1]/cX[i][2];
            Lx[2*i+1][3] = 1+pow(xd[i][1],2);
            Lx[2*i+1][4] = -xd[i][0]*xd[i][1];
            Lx[2*i+1][5] = -xd[i][0];
    }        
}


void tp2DVisualServoingOnePoint()
{

    //-------------------------------------------------------------
    // Mise en oeuvre des courbes
    vpPlot plot(4, 700, 700, 100, 200, "Curves...");

    char title[40];
    strncpy( title, "||e||", 40 );
    plot.setTitle(0,title);
    plot.initGraph(0,1);

    strncpy( title, "x-xd", 40 );
    plot.setTitle(1, title);
    plot.initGraph(1,2);

    strncpy( title, "camera velocity", 40 );
    plot.setTitle(2, title);
    plot.initGraph(2,6);

    strncpy( title, "Camera position", 40 );
    plot.setTitle(3, title);
    plot.initGraph(3,6);

    //-------------------------------------------------------------
    // Affichage des images
    vpImage<unsigned char> I(400,600) ;
    vpDisplayX d ;
    d.init(I) ;
    vpDisplay::display(I);
    vpCameraParameters cam(400,400,300,200) ;

    //-------------------------------------------------------------

    //Definition de la scene
    vpHomogeneousMatrix cTw (0,0,1,  0,0,0) ;

    //position of the point in the world frame
    vpColVector wX(3) ;
    wX[0] = 0.5 ; wX[1] = 0.2 ; wX[2] = -0.5 ;

    // a chaque fois que vous verez size metter la bonne taille à la place
    int size = 2;

    vpColVector e(size) ; //
	e[0] = 10;
	e[1] = 10;

    // position courante, position desiree
    vpColVector x(size), xd(size) ;
    vpColVector cX(4);
    //matrice d'interaction
    vpMatrix Lx(size,6) ;

    // position desirée  (au centre de l'image x=0, y=0)
    xd[0]=0;
    xd[1]=0;

    // vitesse de la camera
    vpColVector v(6) ;
    double lambda = 0.1 ;
    int iter = 0 ;

    while (fabs(e.sumSquare()) > 1e-6)
    {
        // calcul de la position des points dans l'image en fonction de cTw
        changeFrame(wX,cTw,cX);
        
       // instancier x
	    project(cX,x);

        //calcul de l'erreur
		e = x-xd;
        // Calcul de la matrice d'interaction
		computeInteractionMatrix(&cX,&xd,Lx,1);

        //calcul de la loi de commande v= ...

        v = - lambda * Lx.pseudoInverse() * e;

        // Ne pas modifier la suite
        //mise a jour de la position de la camera
        cTw = vpExponentialMap::direct(v).inverse()* cTw ;

        cout << "iter "<< iter <<" : "<< e.t() << endl ;

        iter++ ;

        //mise a jour des courbes
        vpPoseVector ctw(cTw) ;
        plot.plot(0,0,iter, e.sumSquare()) ;
        plot.plot(1,iter, e) ;
        plot.plot(2,iter, v) ;
        plot.plot(3,iter, ctw) ;
        //mise a jour de l'image
        display(cam,I,x,xd) ;
    }

    // sauvegarde des courbes
    plot.saveData(0,"e.txt","#");
    plot.saveData(1,"error.txt","#");
    plot.saveData(2,"v.txt","#");
    plot.saveData(3,"p.txt","#");

    int a ; cin >> a ;
    // sauvegarde de l'image finale
    {
        vpImage<vpRGBa>  Irgb ;
        vpDisplay::getImage(I,Irgb) ;
        vpImageIo::write(Irgb,"1pt.jpg") ;
    }
    cout << "Clicker sur l'image pour terminer" << endl ;
    vpDisplay::getClick(I) ;
}




void tp2DVisualServoingFourPoint()
{
    //-------------------------------------------------------------
    // Mise en oeuvre des courbes
    vpPlot plot(4, 700, 700, 100, 200, "Curves...");

    char title[40];
    strncpy( title, "||e||", 40 );
    plot.setTitle(0,title);
    plot.initGraph(0,1);

    strncpy( title, "x-xd", 40 );
    plot.setTitle(1, title);
    plot.initGraph(1,8);

    strncpy( title, "camera velocity", 40 );
    plot.setTitle(2, title);
    plot.initGraph(2,6);

    strncpy( title, "camera position", 40 );
    plot.setTitle(3, title);
    plot.initGraph(3,6);

    //-------------------------------------------------------------
    // Affichage des images
    vpImage<unsigned char> I(400,600) ;
    vpDisplayX d ;
    d.init(I) ;
    vpDisplay::display(I);
    vpCameraParameters cam(400,400,300,200) ;

    //-------------------------------------------------------------

    //positions initiale (à tester)
    //vpHomogeneousMatrix cTw (-0.2, -0.1, 1.3,
    //                         vpMath::rad(10), vpMath::rad(20), vpMath::rad(30) ) ;
    //vpHomogeneousMatrix cTw (0,0,1.3,  0,0,0) ;
    //vpHomogeneousMatrix cTw (0,0,1,  0,0,vpMath::rad(45)) ;
    //vpHomogeneousMatrix cTw (0, 0, 1,  0, 0, vpMath::rad(90)) ;
    vpHomogeneousMatrix cTw (0, 0, 1,  0, 0, vpMath::rad(180)) ;

    // position finale
    vpHomogeneousMatrix cdTw (0,0,1,  0,0,0) ;

    // position des point dans le repere monde Rw
    vpColVector wX[4] ;
    for (int i = 0 ; i < 4 ; i++) wX[i].resize(3) ;

    double M = 0.3 ;
    wX[0][0] = -M     ; wX[0][1] = -M     ; wX[0][2] = 0 ;
    wX[1][0] = M ;      wX[1][1] = -M     ; wX[1][2] = 0 ;
    wX[2][0] = M ;      wX[2][1] =  M;      wX[2][2] = 0 ;
    wX[3][0] =  -M    ; wX[3][1] =  M;      wX[3][2] = 0 ;

    int size = 4 ;
    vpColVector e(8) ;
    //
    e[0] = 10; e[1] = 10; e[2] = 10; e[3] = 10;
    e[4] = 10; e[5] = 10; e[6] = 10; e[7] = 10;

    vpColVector x[4], xd[4] ;
    vpColVector cX[4];
    for (int i = 0 ; i < 4 ; i++) {
        xd[i].resize(2) ;
        x[i].resize(2);
        cX[i].resize(3);
    }
    vpMatrix Lx(8,6) ;

    //initialisation de la position désire des points dans l'image en fonction de cdTw
    for(int i=0; i<4;i++){
        changeFrame(wX[i],cdTw,cX[i]);
        project(cX[i],xd[i]);
    }

    vpColVector v(6) ;
    double lambda = 0.1 ;
    int iter = 0 ;

    while (fabs(e.sumSquare()) > 1e-6)
    {

        // calcul de la position des points dans l'image en fonction de cTw
        for(int i=0; i<4;i++){
            changeFrame(wX[i],cTw,cX[i]);
            project(cX[i],x[i]);
        }

        // Calcul de la matrice d'interaction
        computeInteractionMatrix(cX,x,Lx,4);

        //calcul de l'erreur
        for (int i = 0; i<4; i++){
            e[2*i] = x[i][0]-xd[i][0];
            e[2*i+1] = x[i][1]-xd[i][1];
        }

        //calcul de la loi de commande
        v = - lambda * Lx.pseudoInverse()*e;

        //mise a jour de la position de la camera
        cTw = vpExponentialMap::direct(v).inverse()* cTw ;

        cout << "iter "<< iter <<" : "<< e.t() << endl ;
        iter++ ;

       //mise a jour des courbes
        vpPoseVector ctw(cTw) ;
        plot.plot(0,0,iter, e.sumSquare()) ;
        plot.plot(1,iter, e) ;
        plot.plot(2,iter, v) ;
        plot.plot(3,iter, ctw) ;
        //mise a jour de l'image
        for(int i = 0;i<4;i++)
        {
            display(cam,I,x[i],xd[i]) ;
        }

    }

    // sauvegarde des courbes
    plot.saveData(0,"e.txt","#");
    plot.saveData(1,"error.txt","#");
    plot.saveData(2,"v.txt","#");
    plot.saveData(3,"p.txt","#");

    // sauvegarde de l'image finale
    {
        vpImage<vpRGBa>  Irgb ;
        vpDisplay::getImage(I,Irgb) ;
        vpImageIo::write(Irgb,"4pt.jpg") ;
    }
    cout << "Clicker sur l'image pour terminer" << endl ;
    vpDisplay::getClick(I) ;
}


void
computeError3D(const vpHomogeneousMatrix &cdTw, const vpHomogeneousMatrix &cTw, vpColVector &error , vpHomogeneousMatrix  &cdTc )
{
    cdTc = cdTw*cTw.inverse();
    vpTranslationVector t = cdTc.getTranslationVector();
    vpThetaUVector thetaU = cdTc.getRotationMatrix().getThetaUVector();
    error[0] = t[0];
    error[1] = t[1];
    error[2] = t[2];
    error[3] = thetaU[0];
    error[4] = thetaU[1];
    error[5] = thetaU[2];
}


void
computeInteractionMatrix3D(const vpHomogeneousMatrix &cdTc, vpMatrix &Lx)
{
    vpMatrix Lw(3,3);
    vpMatrix I(3,3);
    vpRotationMatrix rot = cdTc.getRotationMatrix();
    vpThetaUVector thetaU = rot.getThetaUVector();
    double theta = thetaU.getTheta();
    vpColVector u = thetaU.getU();
    I[0][0] = 1; I[1][1] = 1; I[2][2] = 1;
    Lw = I + theta/2 * vpColVector::skew(u) + (1- vpMath::sinc(theta)/pow(vpMath::sinc(theta/2),2)) * vpColVector::skew(u)*vpColVector::skew(u);
    for (int i = 0; i<6; i++){
        for(int j=0; j<6; j++){
            if(i<3 && j <3){
                Lx[i][j]=rot[i][j];
            }
            else if(i>=3 && j>=3){
                Lx[i][j]=Lw[i-3][j-3];
            }
            else{
                Lx[i][j] = 0;
            }
        }
        
    }
}


void
tp3DVisualServoing()
{

    vpTRACE("begin" ) ;

    vpPlot plot(4, 700, 700, 100, 200, "Curves...");

    char title[40];
    strncpy( title, "||e||", 40 );
    plot.setTitle(0,title);
    plot.initGraph(0,1);

    strncpy( title, "x-xd", 40 );
    plot.setTitle(1, title);
    plot.initGraph(1,6);

    strncpy( title, "camera velocity", 40 );
    plot.setTitle(2, title);
    plot.initGraph(2,6);

    strncpy( title, "Camera position", 40 );
    plot.setTitle(3, title);
    plot.initGraph(3,6);

    //Definition de la scene
    vpHomogeneousMatrix cTw (-0.2, -0.1, 1.3,
                            vpMath::rad(10), vpMath::rad(20), vpMath::rad(30) ) ;
    
    //vpHomogeneousMatrix cTw (0,0,1.3,  0,0,0) ;
    //vpHomogeneousMatrix cTw (0,0,1,  0,0,vpMath::rad(45)) ;
    //vpHomogeneousMatrix cTw (0, 0, 1,  0, 0, vpMath::rad(90)) ;
    //vpHomogeneousMatrix cTw (0, 0, 1,  0, 0, vpMath::rad(180)) ;
    vpHomogeneousMatrix cdTw (0,0,1,  0,0,0) ;

    vpHomogeneousMatrix cdTc;

    int size = 6 ;
    vpColVector e(size) ; //
    e[0] = 10;
    vpMatrix Lx(size,size) ;

    vpColVector v(size) ;
    double lambda = 0.1 ;
    int iter = 0 ;
    while (fabs(e.sumSquare()) > 1e-6)
    {

        // Calcul de l'erreur
        computeError3D(cdTw, cTw, e, cdTc) ;

        // Calcul de la matrice d'interaction
        computeInteractionMatrix3D(cdTc,Lx) ;

        // Calcul de la loi de commande
        v = -lambda * Lx.pseudoInverse()*e;

        // Mis à jour de la position de la camera
        cTw = vpExponentialMap::direct(v).inverse()* cTw ;

        cout << "iter "<< iter <<" : "<< e.t() << endl ;

        iter++ ;

        //mis a jour de courbes
        vpPoseVector crw(cTw) ;
        plot.plot(0,0,iter, e.sumSquare()) ;
        plot.plot(1,iter, e) ;
        plot.plot(2,iter, v) ;
        plot.plot(3,iter,crw) ;
    }

    // sauvegarde des courbes
    plot.saveData(0,"e.txt","#");
    plot.saveData(1,"error.txt","#");
    plot.saveData(2,"v.txt","#");
    plot.saveData(3,"p.txt","#");

    int a ; cin >> a ;

}


void
tp2DVisualServoingFourPointMvt()
{
    //-------------------------------------------------------------
    // Mise en oeuvre des courbes
    vpPlot plot(4, 700, 700, 100, 200, "Curves...");

    char title[40];
    strncpy( title, "||e||", 40 );
    plot.setTitle(0,title);
    plot.initGraph(0,1);

    strncpy( title, "x-xd", 40 );
    plot.setTitle(1, title);
    plot.initGraph(1,8);

    strncpy( title, "camera velocity", 40 );
    plot.setTitle(2, title);
    plot.initGraph(2,6);

    strncpy( title, "camera position", 40 );
    plot.setTitle(3, title);
    plot.initGraph(3,6);

    //-------------------------------------------------------------
    // Affichage des images
    vpImage<unsigned char> I(400,600) ;
    vpDisplayX d ;
    d.init(I) ;
    vpDisplay::display(I);
    vpCameraParameters cam(400,400,300,200) ;

    //-------------------------------------------------------------


    //positions initiale (à tester)
    vpHomogeneousMatrix cTw (-0.2, -0.1, 1.3,
                             vpMath::rad(10), vpMath::rad(20), vpMath::rad(30) ) ;
    //vpHomogeneousMatrix cTw (0,0,1.3,  0,0,0) ;
    //vpHomogeneousMatrix cTw (0,0,1,  0,0,vpMath::rad(45)) ;
    //vpHomogeneousMatrix cTw (0, 0, 1,  0, 0, vpMath::rad(90)) ;
    //vpHomogeneousMatrix cTw (0, 0, 1,  0, 0, vpMath::rad(180)) ;

    // position finale
    vpHomogeneousMatrix cdTw (0,0,1,  0,0,0) ;

    // position des point dans le repere monde Rw
    vpColVector wX[4] ;
    for (int i = 0 ; i < 4 ; i++) wX[i].resize(3) ;

    double M = 0.3 ;
    wX[0][0] = -M     ; wX[0][1] = -M     ; wX[0][2] = 0 ;
    wX[1][0] = M ;      wX[1][1] = -M     ; wX[1][2] = 0 ;
    wX[2][0] = M ;      wX[2][1] =  M;      wX[2][2] = 0 ;
    wX[3][0] =  -M    ; wX[3][1] =  M;      wX[3][2] = 0 ;

    int size = 4 ;
    vpColVector e(8) ;
    //
    e[0] = 10; e[1] = 10; e[2] = 10; e[3] = 10;
    e[4] = 10; e[5] = 10; e[6] = 10; e[7] = 10;

    vpColVector x[4], xd[4] ;
    vpColVector cX[4];
    for (int i = 0 ; i < 4 ; i++) {
        xd[i].resize(2) ;
        x[i].resize(2);
        cX[i].resize(3);
    }
    vpMatrix Lx(8,6) ;

    //initialisation de la position désire des points dans l'image en fonction de cdTw
    for(int i=0; i<4;i++){
        changeFrame(wX[i],cdTw,cX[i]);
        project(cX[i],xd[i]);
    }

    vpColVector v(6) ;
    double lambda = 0.1 ;
    int iter = 0 ;

    vpColVector cumulError(8) ;
    cumulError[0] = 0; cumulError[1] = 0; cumulError[2] = 0; cumulError[3] = 0;
    cumulError[4] = 0; cumulError[5] = 0; cumulError[6] = 0; cumulError[7] = 0;
    double nu = 0.005;

    while (fabs(e.sumSquare()) > 1e-16)
    {
        // les points sont animés d'un mouvement de 1cm/s en x dans Rw
        for (int i = 0 ; i < 4 ; i++) wX[i][0] += 0.01 ;

        // calcul de la position des points dans l'image en fonction de cTw
        for(int i=0; i<4;i++){
            changeFrame(wX[i],cTw,cX[i]);
            project(cX[i],x[i]);
        }

        // Calcul de la matrice d'interaction
        computeInteractionMatrix(cX,xd,Lx,4);

        //calcul de l'erreur
        for (int i = 0; i<4; i++){
            e[2*i] = x[i][0]-xd[i][0];
            e[2*i+1] = x[i][1]-xd[i][1];
        }

        cumulError = cumulError + nu*e;
        //calcul de la loi de commande
        v = Lx.pseudoInverse()*(- lambda * e - cumulError);

        //mise a jour de la position de la camera
        cTw = vpExponentialMap::direct(v).inverse()* cTw ;

        cout << "iter "<< iter <<" : "<< e.t() << endl ;
        iter++ ;

       //mise a jour des courbes
        vpPoseVector ctw(cTw) ;
        plot.plot(0,0,iter, e.sumSquare()) ;
        plot.plot(1,iter, e) ;
        plot.plot(2,iter, v) ;
        plot.plot(3,iter, ctw) ;
        //mise a jour de l'image
        for(int i = 0;i<4;i++)
        {
            display(cam,I,x[i],xd[i]) ;
        }

    }

    // sauvegarde des courbes
    plot.saveData(0,"e.txt","#");
    plot.saveData(1,"error.txt","#");
    plot.saveData(2,"v.txt","#");
    plot.saveData(3,"p.txt","#");

    // sauvegarde de l'image finale
    {
        vpImage<vpRGBa>  Irgb ;
        vpDisplay::getImage(I,Irgb) ;
        vpImageIo::write(Irgb,"4pt.jpg") ;
    }
    cout << "Clicker sur l'image pour terminer" << endl ;
    vpDisplay::getClick(I) ;
}



int main(int argc, char** argv)
{

    //tp2DVisualServoingOnePoint() ;
    //tp2DVisualServoingFourPoint() ;
    //tp3DVisualServoing() ;
    tp2DVisualServoingFourPointMvt() ;

}
