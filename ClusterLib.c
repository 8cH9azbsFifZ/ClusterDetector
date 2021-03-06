/*
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    (C)opyright   Gerolf Ziegenhain 2008 <gerolf@ziegenhain.com>

    Based on the algorithm in S. Stoddard JCP(1978)27:291

*/

#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#include "Numeric/arrayobject.h"


#define EOK 0
#define EMEM 1
#define ESYN 2
#define EFILE 3

#define X 0
#define Y 1
#define Z 2

#define TRUE 1
#define FALSE 0

#define MAGIC_MIN 1e10
#define MAGIC_MAX -1e10

#define MAX(a,b) (a>b)?a:b
#define MIN(a,b) (a>b)?b:a

#define FLOAT float
#define INT int

#define USAGE  "-h     help\n" \
               "-r     rcut\n" \
               "-n     use neighborlist (radius)\n" \
               "-l     use linked lists\n" \
               "-f     infile (name)\n" \
               "-XYZ   peridicity\n"

#define CLK 100000.

#define SQ(a) ((a)*(a))

/* --------------------------------------------------------------------------------------------- */

typedef struct s_linkedlist {
   unsigned int particle;
   struct s_linkedlist *next;
} t_linkedlist;

/* --------------------------------------------------------------------------------------------- */

int Periodic[3];
float Length[3];
unsigned int NCellSide[3], NCell;
float CellLength[3], CellLengthI[3];
t_linkedlist **CellList;
float min[3], max[3];

/* --------------------------------------------------------------------------------------------- */

double SqDist (double pos1, double pos2, int dim) {
   if (Periodic[dim] == TRUE)
      return SQ( fmod(pos1-pos2,Length[dim]) );
   else
      return SQ(pos1-pos2);
}

/* --------------------------------------------------------------------------------------------- */

void DetermineBoxLength (unsigned int n, float *x, float *y, float *z, float *lx, float *ly, float *lz) {
   printf ("Determine boxlength\n");
   min[X] = min[Y] = min[Z] = MAGIC_MIN;
   max[X] = max[Y] = max[Z] = MAGIC_MAX;
   unsigned int i;
   for (i = 0; i < n; i++) {
      min[X] = MIN(min[X],x[i]);
      min[Y] = MIN(min[Y],y[i]);
      min[Z] = MIN(min[Z],z[i]);
      max[X] = MAX(max[X],x[i]);
      max[Y] = MAX(max[Y],y[i]);
      max[Z] = MAX(max[Z],z[i]);
   }

   double off[3];
   off[X] = 0.;
   off[Y] = 0.;
   off[Z] = 0.;
   unsigned int ShiftFlag = FALSE;
   if (min[X] < 0.) {
      ShiftFlag = TRUE;
      off[X] = -min[X];
   }

   if (min[Y] < 0.) {
      ShiftFlag = TRUE;
      off[Y] = -min[Y];
    }
                     
   if (min[Z] < 0.) {
      ShiftFlag = TRUE;
      off[Z] = -min[Z];
   }
                  
   if (ShiftFlag) {
      for (i = 0; i < n; i++) {
         x[i] += off[X];
         y[i] += off[Y];
         z[i] += off[Z];
      }
      unsigned int j;
      for (j = 0; j < 3; j++) {
         min[j] += off[j];
         max[j] += off[j];
      }
   }

   *(lx) = max[X] - min[X];
   *(ly) = max[Y] - min[Y];
   *(lz) = max[Z] - min[Z];
}

/* --------------------------------------------------------------------------------------------- */

void LinkedListInsert (unsigned int particle, t_linkedlist **list) {
   t_linkedlist *new = (t_linkedlist *)calloc (1, sizeof (t_linkedlist));
   if (!new)
      exit (EMEM);
   new->next = *list;
   new->particle = particle;
   *list = new;
}

/* --------------------------------------------------------------------------------------------- */

void LinkedListPrint (t_linkedlist *list) {
   t_linkedlist *ptr = list;
   while (ptr != NULL) {
      ptr = ptr->next;
   }
}

/* --------------------------------------------------------------------------------------------- */

#define CoordinateToCell(x,y,z) (int) (    floor (x*CellLengthI[X]) +    floor (y*CellLengthI[Y])*NCellSide[X] +    floor (z*CellLengthI[Z])*NCellSide[X]*NCellSide[Y]   )


int *BuildLinkedList (int n, float *x, float *y, float *z, float RcNeighbor) {
   printf ("Build linked list\n");
   unsigned int i;

   float CellSize = RcNeighbor; 
   NCellSide[X] = (unsigned int) ceil (Length[X] / CellSize);
   NCellSide[Y] = (unsigned int) ceil (Length[Y] / CellSize);
   NCellSide[Z] = (unsigned int) ceil (Length[Z] / CellSize);
   NCell = NCellSide[X] * NCellSide[Y] * NCellSide[Z];

   CellLength[X] = Length[X] / (float)NCellSide[X];
   CellLength[Y] = Length[Y] / (float)NCellSide[Y];
   CellLength[Z] = Length[Z] / (float)NCellSide[Z];

   CellLengthI[X] = (float)NCellSide[X] / Length[X];
   CellLengthI[Y] = (float)NCellSide[Y] / Length[Y];
   CellLengthI[Z] = (float)NCellSide[Z] / Length[Z];

   CellList = (t_linkedlist **) calloc (NCell, sizeof (t_linkedlist));
   
   for (i = 0; i < n; i++) 
      LinkedListInsert (i, CellList + CoordinateToCell (x[i], y[i], z[i]) );

   return NULL;
}

/* --------------------------------------------------------------------------------------------- */

int *BuildNeighborList (int n, float *x, float *y, float *z, float RcNeighbor) {
   printf ("Build neighbos list\n");
   return BuildLinkedList (n, x, y, z, RcNeighbor);
}

/* --------------------------------------------------------------------------------------------- */

void MoleculeListNeighbor (int Atoms, float *x, float *y, float *z, float RcClusterSq) {
   printf ("Clustering: neighbor list\n");
   double start = clock();
   
   float z0=0.0;
   float RcCluster = sqrt (RcClusterSq);
   int i, j, k, Lk, Size, *List,*Spectrum,s,n;
   
   float xCell[3];
   int m[3];
   t_linkedlist *ptr;

   List = (INT *) malloc (sizeof (INT) * Atoms);
   Spectrum = (INT *) malloc (sizeof (INT) * Atoms+1);
   if (!List || !Spectrum)
      exit (EMEM);

   for (i = 0; i < Atoms; i++)  
      List[i] = i;

   #define NOT_YET_CLUSTERED(i) (i ==  List[i])
   for (i = 0; i < Atoms; i++) {
      if (NOT_YET_CLUSTERED(i)) {
         // set j
         j = i; 
         do {
            // sum over all neighbor cells
            for (m[X] = -1; m[X] <= 1; m[X]++) {
               xCell[X] = x[j]+(float)m[X]*RcCluster;
               if (xCell[X] < min[X] || xCell[X] > max[X]) continue;

               for (m[Y] = -1; m[Y] <= 1; m[Y]++) {
                  xCell[Y] = y[j]+(float)m[Y]*RcCluster;
                  if (xCell[Y] < min[Y] || xCell[Y] > max[Y]) continue;

                  for (m[Z] = -1; m[Z] <= 1; m[Z]++) {
                     xCell[Z] = z[j]+(float)m[Z]*RcCluster;
                     if (xCell[Z] < min[Z] || xCell[Z] > max[Z]) continue;
                     
                     // sum over all particles in cell
                     for (ptr = *(CellList + CoordinateToCell (xCell[X], xCell[Y], xCell[Z])); ptr != NULL; ptr = ptr->next) {
                        k = ptr->particle;
                        //if (k <= i) continue;
                        if (NOT_YET_CLUSTERED(k)) {
                           #define RSQ(j,k) SqDist(x[j],x[k],X) + SqDist(y[j],y[k],Y) + SqDist(z[j],z[k],Z)
                           // FIXME: LinkedList
                           if (RSQ(j,k) <= RcClusterSq) { Lk = List[k]; List[k] = List[j];  List[j] = Lk; }
                        }
                     }

                  }
               }
            }
            // set j
            j = List[j];

         }  while (i != j);
      }  /* if i */
   }  /* for i */
   
   double stop = clock();

   for (i = 0; i < Atoms+1; i++) 
      Spectrum[i] = 0;
   for (i = 0; i < Atoms; i++) {
      if (List[i] >= 0) {
         Size = 1;  j = List[i];  List[i] = -1;
         while (j != i) { Size++;   k = List[j];  List[j] = -1;  j = k; }
         //Bestimmung der Groesse des Clusters im Startteilchen i
         Spectrum[Size] += 1;
      }  /* if List */
   }  /* for i */
   printf("\nCoherent Clusters below %4.4f A\n", RcClusterSq);
   printf("Size Times Atoms\n");
   for (i = 0; i < Atoms+1; i++) // Atoms + 1 because 0 is also a possible Count   
   {
      if ((s = Spectrum[i]) != 0) 
         printf ("%6d %6d %6d\n", i, s, i*s);
   }

   free (List);
   free (Spectrum);
}

/* --------------------------------------------------------------------------------------------- */

void MoleculeList (int Atoms, float *x, float *y, float *z, float RcClusterSq) {
   printf ("Clustering\n");
   double start = clock();
   
   FLOAT     Xj, Yj, Zj, RSq, z0=0.0;

   INT        i, j, k, Lk, Size, *List,*Spectrum,s,n;

   List = (INT *) malloc (sizeof (INT) * Atoms);
   Spectrum = (INT *) malloc (sizeof (INT) * Atoms+1);
   if (!List || !Spectrum)
      exit (EMEM);

   for (i = 0; i < Atoms; i++)  
      List[i] = i;

   for (i = 0; i < Atoms-1; i++) {
      if (i == List[i]) {
         j = i;  Xj = x[j];  Yj = y[j];  Zj = z[j];
         do {
            for (k = i+1; k < Atoms; k++) {
               Lk = List[k];
               if (k == Lk) {
                  RSq = SqDist(Xj,x[k],X) + SqDist(Yj,y[k],Y) + SqDist(Zj,z[k],Z);
                  if (RSq <= RcClusterSq) { List[k] = List[j];  List[j] = Lk; }
               }  /* if k */
            }  /* for k */
            j = List[j];  Xj = x[j];  Yj = y[j];  Zj = z[j];
         }  while (i != j);
      }  /* if i */
   }  /* for i */
   
   double stop = clock();

   for (i = 0; i < Atoms+1; i++) 
      Spectrum[i] = 0;
   for (i = 0; i < Atoms; i++) {
      if (List[i] >= 0) {
         Size = 1;  j = List[i];  List[i] = -1;
         while (j != i) { Size++;   k = List[j];  List[j] = -1;  j = k; }
         //Bestimmung der Groesse des Clusters im Startteilchen i
         Spectrum[Size] += 1;
      }  /* if List */
   }  /* for i */
   printf("\nCoherent Clusters below %4.4f A\n", RcClusterSq);
   printf("Size Times Atoms\n");
   for (i = 0; i < Atoms+1; i++) // Atoms + 1 because 0 is also a possible Count   
   {
      if ((s = Spectrum[i]) != 0) 
         printf ("%6d %6d %6d\n", i, s, i*s);
   }

   free (List);
   free (Spectrum);
}

/* --------------------------------------------------------------------------------------------- */

static PyObject * clusterdetector(PyObject *self, PyObject *args) {
   PyObject *inx, *iny, *inz;
   PyArrayObject *xx, *yy, *zz;
   double rcut;
   int list = 0;

   if (!PyArg_ParseTuple (args, "OOOd|i", &inx, &iny, &inz, &rcut, &list))
      return NULL;

   xx = (PyArrayObject *)PyArray_ContiguousFromObject (inx, PyArray_DOUBLE, 1, 0);
   yy = (PyArrayObject *)PyArray_ContiguousFromObject (iny, PyArray_DOUBLE, 1, 0);
   zz = (PyArrayObject *)PyArray_ContiguousFromObject (inz, PyArray_DOUBLE, 1, 0);

   int n = xx->dimensions[0];
   
   float *x = (float *) calloc (n, sizeof (float));
   float *y = (float *) calloc (n, sizeof (float));
   float *z = (float *) calloc (n, sizeof (float));

   int i;
   for (i=0;i <n;i++) {
      x[i] = (float)*(double*) (xx->data+i*xx->strides[0]);
      y[i] = (float)*(double*) (yy->data+i*yy->strides[0]);
      z[i] = (float)*(double*) (zz->data+i*zz->strides[0]);
   }
   
   DetermineBoxLength (n, x, y, z, &(Length[X]), &(Length[Y]), &(Length[Z]));

   float RcClusterSq, RcCluster = 0.;
   float RcNeighbor = 0;
   RcNeighbor = RcCluster;

   if (list != 0) {
      BuildNeighborList (n, x, y, z, RcNeighbor);
      MoleculeListNeighbor (n, x, y, z, RcClusterSq);
   } else {
      MoleculeList (n, x, y, z, RcClusterSq);
   }


   Py_DECREF (xx);
   Py_DECREF (yy);
   Py_DECREF (zz);

   return  PyFloat_FromDouble (Length[Z]);
}


static struct PyMethodDef cluster_methods[] = {
   {"clusterdetector", clusterdetector, 1},
   {NULL, NULL}
};

PyMODINIT_FUNC initcluster (void) {
   Py_InitModule ("cluster", cluster_methods);
   import_array ();
}
