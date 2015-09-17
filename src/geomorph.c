/* Copyright 2012 Dean C. Adams */

/* This file is part of the R-package `geomorph'. */

/*   C code for Generalized Procrustes Analysis 
    
* This code is designed to perform Generalized Procrustes analysis for points, curves, and surfaces.
* Routines for for OPA and GPA are first, followed by routines for sliding semilandmarks along 
* their tangent directions.
*
*/

#include <R.h>
#include <math.h> 
#include <Rinternals.h>
#include <Rdefines.h>
#include <stdio.h>
#include <R_ext/Lapack.h>
#include <string.h>
#include <float.h>
#include <stdlib.h> 
#include <string.h> 
#include "array.h"

/*  Routines for Generalized Procrustes Analysis */

void OPA(int p2in, int k2in, double RefMat[p2in][k2in], double S2Mat[p2in][k2in], double *OPAres) 
{
  int i, j, p2, k2,kk;  
  p2 = p2in;
  k2 = k2in;
// Regress t(Ref)*S2
  double YtX[k2][k2]; 
   for (i=0;i<k2;i++){
     for (j=0;j<k2;j++){
       YtX[i][j]=0;
     }
   }
   for (i=0;i<k2;i++){
     for (j=0;j<k2;j++){
       for (kk=0;kk<p2;kk++){
         YtX[i][j]=YtX[i][j] +RefMat[kk][i]*S2Mat[kk][j]; //note: ref[kk][i] b/c transpose
       }
     }
   }
// Identify reflections: if det(YtX) < 0, then reflections exist 
   double det;
   if(k2 == 2){det = (YtX[0][0]*YtX[1][1])-(YtX[1][0]*YtX[0][1]);}  
	  else{det = YtX[0][0]*YtX[1][1]*YtX[2][2] + YtX[0][1]*YtX[1][2]*YtX[2][0] +
     YtX[0][2]*YtX[1][0]*YtX[2][1] - YtX[0][0]*YtX[1][2]*YtX[2][1] -
     YtX[0][1]*YtX[1][0]*YtX[2][2] - YtX[0][2]*YtX[1][1]*YtX[2][0];}
    
// SVD of YtX
   double *Umat, *Vtrans, *Dmat; 
   Umat = (double *)R_alloc(k2*k2, sizeof(double));
   Vtrans = (double *)R_alloc(k2*k2, sizeof(double));
   Dmat = (double *)R_alloc(k2, sizeof(double));
   int info,lwork;
   lwork= (6*(k2)*(k2)+6*(k2));
   int *iwork;
   iwork=calloc(8*k2,sizeof(int));     
   double Work[lwork];
   char ch='A';
   F77_CALL(dgesdd)(&ch,&k2,&k2,&YtX[0][0],&k2,Dmat,Vtrans,&k2,Umat,&k2,Work,&lwork,iwork,&info);
  //Assemble H, and project S2*H
   double U[k2][k2],V[k2][k2], H[k2][k2];
   for (i=0;i<k2;i++){
     for (j=0;j<k2;j++){
       U[i][j]=0;
       V[i][j]=0;
       H[i][j]=0;
     }
   }
   for (i=0;i<k2;i++){
     for (j=0;j<k2;j++){
       U[j][i]=Umat[(j*k2)+i];
       V[i][j]=Vtrans[(j*k2)+i];
     }
   }
  // Correct for reflections 
   if(det<0){ 
     for (i=0;i<k2;i++){
	   U[i][k2-1] = -1*U[i][k2-1]; 
     }
   }    
   //Assemble Rotation matrix (H)
   for (i=0;i<k2;i++){
     for (j=0;j<k2;j++){
       for (kk=0;kk<k2;kk++){
         H[i][j]=H[i][j] +V[i][kk]*U[j][kk];  //note: V[j][kk] b/c transpose
       }
     }
   }
   double S3[p2][k2];
   for (i=0;i<p2;i++){
     for (j=0;j<k2;j++){
       S3[i][j]=0;
     }
   }
   for (i=0;i<p2;i++){
     for (j=0;j<k2;j++){
       for (kk=0;kk<k2;kk++){
         S3[i][j]=S3[i][j] +S2Mat[i][kk]*H[kk][j]; 
       }
     }
   }
 memcpy(OPAres,S3,sizeof(double)*p2*k2);

} // end OPA

/*  Main Loop  */
void DoGPA(int *pin, int *kin, int *nin, double *matin, double *res) 
{
  int i,j,jj,mm, p, k,kk,n;  
  p=*pin;
  k=*kin;
  n=*nin;

//Reformat input vector as 3D array
  double ***coords;
  MAKE_3ARRAY(coords,p,k,n);
  double ***coords2;  
  MAKE_3ARRAY(coords2,p,k,n);  
  for (i=0;i<p;i++){
    for (j=0;j<k;j++){
      for (kk=0;kk<n;kk++){
       coords[i][j][kk]=matin[(kk*p*k)+(j*p)+i];
       coords2[i][j][kk]=coords[i][j][kk]; 
     }
    }
  }
//Translate specimens to origin
  double mncoord[k];
  for (i=0;i<n;i++){
    for (kk=0;kk<k;kk++){
      mncoord[kk]=0;
    }
    for (j=0;j<p;j++){
      for (kk=0;kk<k;kk++){
        mncoord[kk]=mncoord[kk]+ coords[j][kk][i];
      }
    }
    for (kk=0;kk<k;kk++){
      mncoord[kk]=mncoord[kk]/p;
    }
    for (j=0;j<p;j++){
      for (kk=0;kk<k;kk++){
        coords[j][kk][i]=coords[j][kk][i]-mncoord[kk];
      }
    }
  }
//Scale specimens to unit centroid-size
 double csize;
  for (i=0;i<n;i++){
    csize=0;
    for (j=0;j<p;j++){
      for (kk=0;kk<k;kk++){
        csize=csize+ (coords[j][kk][i]*coords[j][kk][i]);
      }
    }
    csize=sqrt(csize);
    for (j=0;j<p;j++){
      for (kk=0;kk<k;kk++){
        coords[j][kk][i]=coords[j][kk][i]/csize;
      }
    }
  }
// GPA  setup
  int iter, maxiter;
  double Q,Q1,Q2,minChange;
  double ref[p][k],S2[p][k]; 
  double **DistQ;
  MAKE_2ARRAY(DistQ,n,n);

  minChange = 0.0001;
  Q=0; Q1=0; Q2=0;
  iter=0;  maxiter=5;
  double *OPAout;
  OPAout = (double *)R_alloc(p*k, sizeof(double));
//Initial Q
  for (j=0;j<n;j++){
    for (kk=0;kk<n;kk++){
      DistQ[j][kk]=0;
    }
  }
  for (j=0;j<n;j++){
    for (kk=0;kk<n;kk++){
      for (jj=0;jj<p;jj++){
        for (mm=0;mm<k;mm++){
      DistQ[j][kk]=DistQ[j][kk]+ (coords[jj][mm][j]-coords[jj][mm][kk])*(coords[jj][mm][j]-coords[jj][mm][kk]);
        }
      }
    }
  }
  for (j=0;j<n;j++){
    for (kk=0;kk<n;kk++){
      Q1=Q1+sqrt(DistQ[j][kk]);
    }
  }
  Q1=Q1/2;
  Q=Q1;
// Initial Rotation to Spec 1=Ref
  for (j=0;j<p;j++){
    for (kk=0;kk<k;kk++){
      ref[j][kk]=0;
    }
  }
  for (j=0;j<p;j++){
    for (kk=0;kk<k;kk++){
      ref[j][kk]=coords[j][kk][0]; 
    }
  }
  for (i=0;i<n;i++){
    for (j=0;j<p;j++){
      for (kk=0;kk<k;kk++){
        S2[j][kk]=0;
      }
    }
    for (j=0;j<p;j++){
      for (kk=0;kk<k;kk++){
        S2[j][kk]=coords[j][kk][i]; 
      }
    } 
    for (j=0;j<p*k;j++){
      OPAout[j]=0;
    }
    OPA(p, k, ref, S2, OPAout); 
    for (j=0;j<p;j++){
      for (kk=0;kk<k;kk++){
        coords[j][kk][i]=OPAout[(j*k)+kk];          
      }
    }
  }
  for (j=0;j<n;j++){
    for (kk=0;kk<n;kk++){
      DistQ[j][kk]=0;
    }
  }
  for (j=0;j<n;j++){
    for (kk=0;kk<n;kk++){
      for (jj=0;jj<p;jj++){
        for (mm=0;mm<k;mm++){
      DistQ[j][kk]=DistQ[j][kk]+ (coords[jj][mm][j]-coords[jj][mm][kk])*(coords[jj][mm][j]-coords[jj][mm][kk]);
        }
      }
    }
  }
  for (j=0;j<n;j++){
    for (kk=0;kk<n;kk++){
      Q2=Q2+sqrt(DistQ[j][kk]);
    }
  }
  Q2=Q2/2;
  if(Q1 < Q2){Q=Q1-Q2;}  
  //new check for 1st iteration. If worse, dump and keep old. 9/2015
  if (Q2 > Q1){
    for (i=0;i<p;i++){
      for (j=0;j<k;j++){
        for (kk=0;kk<n;kk++){
          coords[i][j][kk]=coords2[i][j][kk];  
        }
      }
    }
    Q=0;}

// Rotation iterations
while (Q>minChange){
  for (j=0;j<p;j++){
    for (kk=0;kk<k;kk++){
      ref[j][kk]=0;
    }
  }
 for (i=0;i<n;i++){ 
    for (j=0;j<p;j++){
      for (kk=0;kk<k;kk++){
        ref[j][kk]=ref[j][kk]+coords[j][kk][i]/n; 
      }
    }
  }
  for (i=0;i<n;i++){
    for (j=0;j<p;j++){
      for (kk=0;kk<k;kk++){
        S2[j][kk]=0;
      }
    }
    for (j=0;j<p;j++){
      for (kk=0;kk<k;kk++){
        S2[j][kk]=coords[j][kk][i]; 
      }
    } 
    for (j=0;j<p*k;j++){
      OPAout[j]=0;
    }
    OPA(p, k, ref, S2, OPAout); 
    for (j=0;j<p;j++){
      for (kk=0;kk<k;kk++){
        coords[j][kk][i]=OPAout[(j*k)+kk];          
      }
    }
  }
  for (j=0;j<n;j++){
    for (kk=0;kk<n;kk++){
      DistQ[j][kk]=0;
    }
  }
  for (j=0;j<n;j++){
    for (kk=0;kk<n;kk++){
      for (jj=0;jj<p;jj++){
        for (mm=0;mm<k;mm++){
      DistQ[j][kk]=DistQ[j][kk]+ (coords[jj][mm][j]-coords[jj][mm][kk])*(coords[jj][mm][j]-coords[jj][mm][kk]);
        }
      }
    }
  }
  Q2=0;
  for (j=0;j<n;j++){
   for (kk=0;kk<n;kk++){
      Q2=Q2+sqrt(DistQ[j][kk]);
    }
  }
  Q2=Q2/2;
  Q=Q1-Q2;
  Q1=Q2;
  iter=iter+1;
  if(iter==maxiter){Q=0;}
} //End while loop
 //Re-format results into *res
 for (i=0;i<n;i++){
    for (j=0;j<p;j++){
      for (kk=0;kk<k;kk++){
        res[(i*p*k)+(j*k)+kk]=coords[j][kk][i];
      }
    }
  }  
  FREE_3ARRAY(coords);
  FREE_3ARRAY(coords2);  
  FREE_2ARRAY(DistQ);
} //end GPA


/*  GPA 1 Iteration Function  */
void DoGPA1(int *pin, int *kin, int *nin, double *matin, double *res) 
{
  int i,j,p, k,kk,n;  
  p=*pin;
  k=*kin;
  n=*nin;
//Reformat input vector as 3D array
  double ***coords;
  MAKE_3ARRAY(coords,p,k,n);
  for (i=0;i<p;i++){
    for (j=0;j<k;j++){
      for (kk=0;kk<n;kk++){
       coords[i][j][kk]=matin[(kk*p*k)+(j*p)+i];
     }
    }
  }
//Translate specimens to origin
  double mncoord[k];
  for (i=0;i<n;i++){
    for (kk=0;kk<k;kk++){
      mncoord[kk]=0;
    }
    for (j=0;j<p;j++){
      for (kk=0;kk<k;kk++){
        mncoord[kk]=mncoord[kk]+ coords[j][kk][i];
      }
    }
    for (kk=0;kk<k;kk++){
      mncoord[kk]=mncoord[kk]/p;
    }
    for (j=0;j<p;j++){
      for (kk=0;kk<k;kk++){
        coords[j][kk][i]=coords[j][kk][i]-mncoord[kk];
      }
    }
  }
//Scale specimens to unit centroid-size
 double csize;
  for (i=0;i<n;i++){
    csize=0;
    for (j=0;j<p;j++){
      for (kk=0;kk<k;kk++){
        csize=csize+ (coords[j][kk][i]*coords[j][kk][i]);
      }
    }
    csize=sqrt(csize);
    for (j=0;j<p;j++){
      for (kk=0;kk<k;kk++){
        coords[j][kk][i]=coords[j][kk][i]/csize;
      }
    }
  }
// GPA  setup
  double ref[p][k],S2[p][k];
  double *OPAout;
  OPAout = (double *)R_alloc(p*k, sizeof(double));
// Initial Rotation to Spec 1=Ref
  for (j=0;j<p;j++){
    for (kk=0;kk<k;kk++){
      ref[j][kk]=0;
    }
  }
  for (j=0;j<p;j++){
    for (kk=0;kk<k;kk++){
      ref[j][kk]=coords[j][kk][0]; 
    }
  }
  for (i=0;i<n;i++){
    for (j=0;j<p;j++){
      for (kk=0;kk<k;kk++){
        S2[j][kk]=0;
      }
    }
    for (j=0;j<p;j++){
      for (kk=0;kk<k;kk++){
        S2[j][kk]=coords[j][kk][i]; 
      }
    } 
    for (j=0;j<p*k;j++){
      OPAout[j]=0;
    }
    OPA(p, k, ref, S2, OPAout); 
    for (j=0;j<p;j++){
      for (kk=0;kk<k;kk++){
        coords[j][kk][i]=OPAout[(j*k)+kk];          
      }
    }
  }
 //Re-format results into *res
 for (i=0;i<n;i++){
    for (j=0;j<p;j++){
      for (kk=0;kk<k;kk++){
        res[(i*p*k)+(j*k)+kk]=coords[j][kk][i];
      }
    }
  }     
  FREE_3ARRAY(coords);
} //end GPA1

/*  Sliding Semilandmarks for GPA */

int compare_doubles (const void *X, const void *Y){
  double * x = *((double **)X);
  double * y = *((double **)Y);
  if (*x > *y){return 1;}
  else{if (*x < *y){return -1;}
  else{return 0;}
  }
}

void DoSlide(int *ProcDin, int *pin, int *kin, int *nin, double *matin, double *refin, double *res, int *ncurvein, int *curvein, 
	int *nsurfin, int *surfin) 
{
  int i,j,p,ii, k,kk,m,jj,n,ProcD,iter;
  p=*pin;
  k=*kin;
  n=*nin;
  ProcD=*ProcDin;  // 1 = ProcD sliding, 0 = BE sliding
  iter=3;
  double **U;
  double ***coords;
  MAKE_3ARRAY(coords,p,k,n);
  for (i=0;i<p;i++){
    for (j=0;j<k;j++){
      for (kk=0;kk<n;kk++){
       coords[i][j][kk]=matin[(kk*p*k)+(j*p)+i];
     }
    }
  }
  double ref[p][k];
  for (i=0;i<p;i++){
    for (j=0;j<k;j++){
     ref[i][j]=refin[(j*p)+i];
    }
  }
  int ncurve;
  ncurve = *ncurvein;
  int **curves=0;
  if (ncurve >0){
  MAKE_2ARRAY(curves,ncurve,3);
  for (i=0;i<3;i++){
      for (j=0;j<ncurve;j++){
        curves[j][i]=curvein[(i*ncurve)+j];
      }
    }
  } 
  int nsurf;
  nsurf= *nsurfin;
  int *surf=0;
  if (nsurf >0){
  MAKE_1ARRAY(surf,nsurf);
  for (i=0;i<nsurf;i++){
      surf[i]=surfin[i];
    }
  }
  int ncolU;  // Initialize U  
  ncolU=p;
  if(nsurf>0){ ncolU= 2*p; }  
  MAKE_2ARRAY(U,(p*k),ncolU);
  for (ii=0;ii<iter;ii++){
    for (i=0;i<n;i++){
	  for (j=0;j<p*k;j++){
        for (kk=0;kk<ncolU;kk++){
          U[j][kk]=0; 
        }
      }		   
 	  if(ncurve > 0){	   //Curve portion of U
	    for (j=0;j<ncurve;j++){
          U[(curves[j][1]-1)][(curves[j][1]-1)] = coords[(curves[j][2]-1)][0][i]- coords[(curves[j][0]-1)][0][i];
          U[(p+curves[j][1]-1)][(curves[j][1]-1)] = coords[(curves[j][2]-1)][1][i]- coords[(curves[j][0]-1)][1][i];
          if(k==3){
            U[(2*p+curves[j][1]-1)][(curves[j][1]-1)] = coords[(curves[j][2]-1)][2][i]- coords[(curves[j][0]-1)][2][i];
	      }    
        }
      }  
      double NNcoords[5][k], NNmn[k];  
      double *distvec, **distmat;
      MAKE_2ARRAY(distmat,p,p);
      MAKE_1ARRAY(distvec,p);
      double * distvec_p[p];
      if(nsurf >0){ //Surface portion of U: find nearest neigbors (from dist.mat), PCA for tangent directions
        for (j=0;j<p;j++){
          for (kk=0;kk<p;kk++){
            distmat[j][kk]=0;     
          }
        }
	    for (j=0;j<p;j++){
          for (kk=0;kk<p;kk++){
	        for (m=0;m<k;m++){  
              distmat[j][kk]=distmat[j][kk]+((coords[j][m][i]-coords[kk][m][i])*(coords[j][m][i]-coords[kk][m][i]));   
            }
          }
        }
        for (j=0;j<p;j++){
          for (kk=0;kk<p;kk++){
            distmat[j][kk]=sqrt(distmat[j][kk]);     
          }
        }      
	    for (j=0;j<nsurf;j++){
	      for (kk=0;kk<p;kk++){
            distvec[kk]=0;
		  } 
		  for (kk=0;kk<p;kk++){
            distvec[kk]=distmat[(surf[j]-1)][kk];
            distvec_p[kk]=&distvec[kk];
		  } 
          qsort(&distvec_p, p, sizeof(double*), compare_doubles);    //sort distance to find NN. Put coords in matrix, and column-center 
		  double intvec[p];
	      for(kk=0;kk<p;kk++){
		    intvec[kk] = (double)(distvec_p[kk] - &distvec[0]);
	      }
          for (kk=0;kk<5;kk++){
            for (m=0;m<k;m++){
              NNcoords[kk][m]=0;     
            }
          }
          for (kk=0;kk<5;kk++){
            for (m=0;m<k;m++){
              NNcoords[kk][m]=coords[(int)intvec[kk]][m][i];     
            }
          }
          for (m=0;m<k;m++){
            NNmn[m]=0;     
          }
          for (kk=0;kk<5;kk++){
            for (m=0;m<k;m++){
              NNmn[m]=NNmn[m]+NNcoords[kk][m];     
            }
          }
          for (m=0;m<k;m++){
            NNmn[m]=NNmn[m]/5;     
          }
          for (kk=0;kk<5;kk++){
            for (m=0;m<k;m++){
              NNcoords[kk][m]=NNcoords[kk][m]-NNmn[m];     
            }
          }        
          double vcv[k][k];  //SVD on VCV matrix from NN
          for (kk=0;kk<k;kk++){
             for (m=0;m<k;m++){
               vcv[kk][m]=0;     
             }
          }         
          for (jj=0;jj<k;jj++){
            for (kk=0;kk<k;kk++){
              for (m=0;m<5;m++){
                vcv[jj][kk]=vcv[jj][kk]+ NNcoords[m][jj]*NNcoords[m][kk]; 
              }
            }
          }        
          for (kk=0;kk<k;kk++){
            for (m=0;m<k;m++){
              vcv[kk][m]=vcv[kk][m]/4;     
            }
          }         
          double *Umat1, *Vtrans1, *Dmat1; 
          Umat1 = (double *)R_alloc(k*k, sizeof(double));
          Vtrans1 = (double *)R_alloc(k*k, sizeof(double));
          Dmat1 = (double *)R_alloc(k, sizeof(double));
          int info1,lwork1;
          lwork1= (6*(k)*(k)+6*(k));
          int *iwork1;
          iwork1=calloc(8*k,sizeof(int));   
          double Work1[lwork1];
          char ch='A';
          F77_CALL(dgesdd)(&ch,&k,&k,&vcv[0][0],&k,Dmat1,Umat1,&k,Vtrans1,&k,Work1,&lwork1,iwork1,&info1);
          double PCvec[k][k];  //Assemble PCvectors, dump into U
          for (jj=0;jj<k;jj++){
            for (kk=0;kk<k;kk++){
              PCvec[jj][kk]=0;
            }
          }
          for (jj=0;jj<k;jj++){
            for (kk=0;kk<k;kk++){
              PCvec[kk][jj]=Vtrans1[(kk*k)+jj];
            }
          }
          U[surf[j]-1][surf[j]-1] = PCvec[0][0];
          U[(p+surf[j]-1)][surf[j]-1] = PCvec[1][0];
          if(k==3){
            U[(2*p+surf[j]-1)][surf[j]-1] = PCvec[2][0];
          }
          U[surf[j]-1][p+surf[j]-1] = PCvec[0][1];
          U[(p+surf[j]-1)][p+surf[j]-1] = PCvec[1][1];
          if(k==3){
            U[(2*p+surf[j]-1)][p+surf[j]-1] = PCvec[2][1];
          }
        }    // end j loop
      }  //end if surfin loop
      double **Usq, *sumUsq;   //standardize U
      MAKE_2ARRAY(Usq,(p*k),ncolU);
      MAKE_1ARRAY(sumUsq,ncolU);     
      for (j=0;j<(p*k);j++){
        for (kk=0;kk<ncolU;kk++){
          Usq[j][kk]=0; 
        }
      }      
      for (j=0;j<(p*k);j++){
        for (kk=0;kk<ncolU;kk++){
          Usq[j][kk]=U[j][kk]*U[j][kk]; 
        }
      }      
      for (j=0;j<ncolU;j++){
        sumUsq[j]=0; 
      } 
      for (j=0;j<ncolU;j++){
        for (kk=0;kk<(p*k);kk++){
          sumUsq[j]= sumUsq[j]+Usq[kk][j]; 
        }
      }  
      for (j=0;j<ncolU;j++){
        sumUsq[j]=sqrt(sumUsq[j]); 
      } 
      for (j=0;j<(p*k);j++){
        for (kk=0;kk<ncolU;kk++){
	     if (sumUsq[kk] ==0){U[j][kk]=0;}
           else{U[j][kk]=U[j][kk]/sumUsq[kk]; }	        
        }
      }      
      FREE_2ARRAY(Usq);
      FREE_1ARRAY(sumUsq);            
      double *tmpvec, *tmpdiff;
      MAKE_1ARRAY(tmpvec,p*k);
      MAKE_1ARRAY(tmpdiff,p*k);  
      for (j=0;j<(p*k);j++){
        tmpvec[j]=0;
        tmpdiff[j]=0;
      }
      for (j=0;j<p;j++){
        for (kk=0;kk<k;kk++){
         tmpvec[kk*p+j] = coords[j][kk][i]; 
         tmpdiff[kk*p+j] = coords[j][kk][i]- ref[j][kk]; 
        }
      }      
      double **Tpart, *slidvec;
      MAKE_2ARRAY(Tpart,(p*k),(p*k));
      MAKE_1ARRAY(slidvec,p*k);   
     if(ProcD==1){  //ProcD method: Tpart = U%*%t(U)
       for (j=0;j<(p*k);j++){
         for (kk=0;kk<(p*k);kk++){
           Tpart[j][kk]=0; 
         }
       } 
       for (jj=0;jj<(p*k);jj++){
         for (kk=0;kk<(p*k);kk++){
           for (m=0;m<ncolU;m++){   
             Tpart[jj][kk]=Tpart[jj][kk] + U[jj][m]*U[kk][m]; 
           }
         }
       }    
     }   
     double L[p+k+1][p+k+1];
     if(ProcD==0){  //Bending Energy Tpart = U%*% (ginv(t(U)%*%L.inv%*%U) %*%t(U) %*%L.inv)
       double **Pdist;
       MAKE_2ARRAY(Pdist,p,p);
       for (j=0;j<(p*k);j++){
         for (kk=0;kk<(p*k);kk++){
           Tpart[j][kk]=0; 
         }
       } 
       for (j=0;j<(p);j++){
         for (kk=0;kk<(p);kk++){
           Pdist[j][kk]=0; 
         }
       }
       for (j=0;j<(p+k+1);j++){
         for (kk=0;kk<(p+k+1);kk++){
           L[j][kk]=0; 
         }
       }          
  	   for (j=0;j<p;j++){
         for (kk=0;kk<p;kk++){
	       for (m=0;m<k;m++){  
             Pdist[j][kk]=Pdist[j][kk]+((ref[j][m]-ref[kk][m])*(ref[j][m]-ref[kk][m]));   
           }
         }
       }
       for (j=0;j<p;j++){
         for (kk=0;kk<p;kk++){
           Pdist[j][kk]=sqrt(Pdist[j][kk]);     
         }
       }     
       for (j=0;j<p;j++){
         for (kk=0;kk<p;kk++){
	       if(k==3){L[j][kk]=Pdist[j][kk]; }
	       else{ L[j][kk]=(Pdist[j][kk]*Pdist[j][kk])*log(Pdist[j][kk]*Pdist[j][kk]);}  
         }
       } 
       FREE_2ARRAY(Pdist);       
       for (j=0;j<p;j++){
         L[j][j]=0;  // force diag = 0.0 for 2D data
       } 
       for (j=0;j<p;j++){
         L[j][p]=1;  
         L[p][j]=1;    
       }        
       for (j=0;j<p;j++){  
         for (kk=0;kk<k;kk++){
           L[j][p+1+kk]=ref[j][kk];  
           L[p+1+kk][j]=ref[j][kk];    
         }
       }  
       double *Umat2, *Vtrans2, *Dmat2; 
       Umat2 = (double *)R_alloc((p+k+1)*(p+k+1), sizeof(double));
       Vtrans2 = (double *)R_alloc((p+k+1)*(p+k+1), sizeof(double));
       Dmat2 = (double *)R_alloc((p+k+1), sizeof(double));
       int info2,lwork2;
       lwork2=(6*(p+k+1)*(p+k+1)+6*(p+k+1));
       int *iwork2, Lsize;
       Lsize=(p+k+1);
       iwork2=calloc(8*(p+k+1),sizeof(int));   
       double Work2[lwork2];
       char ch1='A';
       F77_CALL(dgesdd)(&ch1,&Lsize,&Lsize,&L[0][0],&Lsize,Dmat2,Umat2,&Lsize,Vtrans2,&Lsize,Work2,&lwork2,iwork2,&info2);
       double **U2, **V2, **Lsolve, **Ltemp, **D2, **BigLinv;
       MAKE_2ARRAY(U2,(p+k+1),(p+k+1));
       MAKE_2ARRAY(V2,(p+k+1),(p+k+1));
       MAKE_2ARRAY(Lsolve,(p+k+1),(p+k+1));
       MAKE_2ARRAY(Ltemp,(p+k+1),(p+k+1));
       MAKE_2ARRAY(D2,(p+k+1),(p+k+1));                           
       MAKE_2ARRAY(BigLinv,(p*k),(p*k));
       for (j=0;j<(p+k+1);j++){
         for (kk=0;kk<(p+k+1);kk++){
           U2[j][kk]=0;
           V2[j][kk]=0;
           D2[j][kk]=0;
           Lsolve[j][kk]=0;
           Ltemp[j][kk]=0;
         }
       }
       for (j=0;j<(p+k+1);j++){
         for (kk=0;kk<(p+k+1);kk++){
           U2[kk][j]=Umat2[(kk*(p+k+1))+j];
           V2[j][kk]=Vtrans2[(kk*(p+k+1))+j];
         }
       }
       for (j=0;j<(p+k+1);j++){
	     D2[j][j]=1/Dmat2[j];  
       }
       for (j=0;j<(p+k+1);j++){
         for (kk=0;kk<(p+k+1);kk++){
           for (m=0;m<(p+k+1);m++){
             Ltemp[kk][j]=Ltemp[kk][j]+U2[m][j]*D2[kk][m];
           }
         }
       } 
       for (j=0;j<(p+k+1);j++){
         for (kk=0;kk<(p+k+1);kk++){
           for (m=0;m<(p+k+1);m++){
             Lsolve[kk][j]=Lsolve[kk][j]+Ltemp[m][j]*V2[m][kk];  //R & C reversed for V2 b/c transpose
           }   
         }
       } 
       for (j=0;j<(p*k);j++){
         for (kk=0;kk<(p*k);kk++){
           BigLinv[j][kk]=0; 
         }
       }
       for (j=0;j<p;j++){
         for (kk=0;kk<p;kk++){
           BigLinv[j][kk]=Lsolve[j][kk];
           BigLinv[j+p][kk+p]=Lsolve[j][kk]; 
           if(k==3){ BigLinv[j+2*p][kk+2*p]=Lsolve[j][kk]; }
         }
       }
       FREE_2ARRAY(U2);
       FREE_2ARRAY(V2);
       FREE_2ARRAY(Lsolve);
       FREE_2ARRAY(Ltemp);
       FREE_2ARRAY(D2);            
       double ULinvU[ncolU][ncolU];       
       double **ULinv; //**ULinvU;   // Tpart for BE = U%*% (ginv(t(U)%*%L.inv%*%U) %*%t(U) %*%L.inv
       MAKE_2ARRAY(ULinv,ncolU,(p*k));
       for (j=0;j<ncolU;j++){
         for (kk=0;kk<(p*k);kk++){
           ULinv[j][kk]=0; 
         }
       }           
       for (j=0;j<ncolU;j++){
         for (kk=0;kk<ncolU;kk++){
           ULinvU[j][kk]=0; 
         }
       }      
       for (j=0;j<ncolU;j++){  //assemble matrix for 'ginv'
         for (kk=0;kk<(p*k);kk++){
           for (m=0;m<(p*k);m++){
             ULinv[j][kk]=ULinv[j][kk]+U[m][j]*BigLinv[kk][m];
           }
         }
       } 
       for (j=0;j<ncolU;j++){  
         for (kk=0;kk<ncolU;kk++){
           for (m=0;m<(p*k);m++){
             ULinvU[j][kk]=ULinvU[j][kk]+ULinv[j][m]*U[m][kk];
           }
         }
       } 
       FREE_2ARRAY(ULinv);    
       double *Umat3, *Vtrans3, *Dmat3; 
       Umat3 = (double *)R_alloc((ncolU)*(ncolU), sizeof(double));
       Vtrans3 = (double *)R_alloc((ncolU)*(ncolU), sizeof(double));
       Dmat3 = (double *)R_alloc((ncolU), sizeof(double));
       int info3,lwork3;
       lwork3=(6*(ncolU)*(ncolU)+6*(ncolU));
       int *iwork3, Lsize3;
       Lsize3=(ncolU);
       iwork3=calloc(8*(ncolU),sizeof(int));   
       double Work3[lwork3];
       char ch3='A';
       F77_CALL(dgesdd)(&ch3,&Lsize3,&Lsize3,&ULinvU[0][0],&Lsize3,Dmat3,Umat3,&Lsize3,Vtrans3,&Lsize3,Work3,&lwork3,iwork3,&info3);
       double tol;
       double **U3,**V3,**ginv, **ginvtemp,**D3;
       MAKE_2ARRAY(U3,ncolU,ncolU);
       MAKE_2ARRAY(V3,ncolU,ncolU);
       MAKE_2ARRAY(D3,ncolU,ncolU);       
       MAKE_2ARRAY(ginv,ncolU,ncolU);
       MAKE_2ARRAY(ginvtemp,ncolU,ncolU);       
       tol= 0.000000001;
       for (j=0;j<(ncolU);j++){
         for (kk=0;kk<(ncolU);kk++){
           U3[j][kk]=0;
           V3[j][kk]=0;
           D3[j][kk]=0;
           ginv[j][kk]=0;
           ginvtemp[j][kk]=0;
         }
       }
       for (j=0;j<(ncolU);j++){
         for (kk=0;kk<(ncolU);kk++){
           U3[kk][j]=Umat3[(kk*(ncolU))+j];
           V3[j][kk]=Vtrans3[(kk*(ncolU))+j];
         }
       }
       for (j=0;j<(ncolU);j++){
	     if (Dmat3[j] > tol){D3[j][j]=1/Dmat3[j];}
           else{D3[j][j]=0;}
       }
       for (j=0;j<(ncolU);j++){
         for (kk=0;kk<(ncolU);kk++){
           for (m=0;m<(ncolU);m++){
             ginvtemp[kk][j]=ginvtemp[kk][j]+U3[m][j]*D3[kk][m];
           }
         }
       } 
       for (j=0;j<(ncolU);j++){
         for (kk=0;kk<(ncolU);kk++){
           for (m=0;m<(ncolU);m++){
             ginv[kk][j]=ginv[kk][j]+ginvtemp[m][j]*V3[m][kk];  //R & C reversed for V3 b/c transpose
           }   
         }
       } 
       FREE_2ARRAY(U3);
       FREE_2ARRAY(V3);
       FREE_2ARRAY(D3);       
       FREE_2ARRAY(ginvtemp); 
       double **ginvU,**ginvULinv;
       MAKE_2ARRAY(ginvU,ncolU,(p*k));
       MAKE_2ARRAY(ginvULinv,ncolU,(p*k));
       for (j=0;j<ncolU;j++){
         for (kk=0;kk<(p*k);kk++){
           ginvU[j][kk]=0; 
           ginvULinv[j][kk]=0;
         }
       } 
       for (j=0;j<ncolU;j++){  
         for (kk=0;kk<(p*k);kk++){
           for (m=0;m<ncolU;m++){
             ginvU[j][kk]=ginvU[j][kk]+ginv[j][m]*U[kk][m];
           }
         }
       } 
       for (j=0;j<ncolU;j++){    
         for (kk=0;kk<(p*k);kk++){
           for (m=0;m<(p*k);m++){
             ginvULinv[j][kk]=ginvULinv[j][kk]+ginvU[j][m]*BigLinv[m][kk];
           }
         }
       }       
       for (j=0;j<(p*k);j++){    
         for (kk=0;kk<(p*k);kk++){
           for (m=0;m<ncolU;m++){
             Tpart[j][kk]=Tpart[j][kk]+U[j][m]*ginvULinv[m][kk];
           }
         }
       }              
       for (j=0;j<(p*k);j++){    
         for (kk=0;kk<(p*k);kk++){
           Tpart[j][kk]=Tpart[j][kk]/2;  // re-scale 
         }
       }
     FREE_2ARRAY(ginv);
     FREE_2ARRAY(BigLinv);   
     FREE_2ARRAY(ginvU);
     FREE_2ARRAY(ginvULinv);          
     }  //ProcD=0 loop         
      for (j=0;j<(p*k);j++){ //Do sliding   
         slidvec[j]=0;
       }           
       for (jj=0;jj<(p*k);jj++){
         for (kk=0;kk<(p*k);kk++){
           slidvec[jj]= slidvec[jj]+ Tpart[jj][kk]*tmpdiff[kk]; 
         }
       }      
       for (jj=0;jj<(p*k);jj++){
         tmpvec[jj]= tmpvec[jj]-slidvec[jj]; 
       }      
       for (j=0;j<k;j++){
         for (kk=0;kk<p;kk++){
           coords[kk][j][i] = tmpvec[kk+p*j]; 
         }
       }     
     FREE_2ARRAY(distmat);
     FREE_1ARRAY(distvec);
     FREE_2ARRAY(Tpart);
     FREE_1ARRAY(slidvec);   
     FREE_1ARRAY(tmpvec);
     FREE_1ARRAY(tmpdiff);  
   }// end i [n] loop
   double tmpref[p][k];  //find new ref after sliding
    for (j=0;j<p;j++){
      for (kk=0;kk<k;kk++){
        tmpref[j][kk]=0;
      }
    }
    for (i=0;i<n;i++){
      for (j=0;j<p;j++){
        for (kk=0;kk<k;kk++){
          tmpref[j][kk]=tmpref[j][kk]+ coords[j][kk][i];
        }
      }
    }  
    for (j=0;j<p;j++){
      for (kk=0;kk<k;kk++){
        tmpref[j][kk]=tmpref[j][kk]/n;
      }
    } 
  }  //end ii loop 
  for (i=0;i<n;i++){  //Re-format results into *res
    for (j=0;j<p;j++){
      for (kk=0;kk<k;kk++){
        res[(i*p*k)+(j*k)+kk]=coords[j][kk][i];
      }
    }
  }
  FREE_2ARRAY(U);
  FREE_3ARRAY(coords);
  FREE_2ARRAY(curves);
  FREE_1ARRAY(surf);
} // end DoSlide    
