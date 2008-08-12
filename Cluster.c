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
                  Christian Anders 2008 <anders@physik.uni-kl.de>

    Based on the algorithm in S. Stoddard JCP(1978)27:291

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#include <string.h>

#include "H5LT.h"
#include "H5TA.h"
#include <hdf5.h>
#undef MAX
#undef MIN


int getopt(int argc, char * const argv[], const char *optstring);
extern char *optarg;
extern int optind, opterr, optopt;
FILE *popen(const char *command, const char *type);


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
#define PROGRESS(n) if (n % 1000 == 0) llog (" %d", n)

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

void llog (const char *format, ...) {
   va_list arg;
   va_start (arg, format);
   vfprintf (stderr, format, arg);
   fflush (stderr);
   va_end (arg);
}

/* --------------------------------------------------------------------------------------------- */

double SqDist (double pos1, double pos2, int dim) {
   if (Periodic[dim] == TRUE)
      return SQ( fmod(pos1-pos2,Length[dim]) );
   else
      return SQ(pos1-pos2);
}

/* --------------------------------------------------------------------------------------------- */

void DetermineBoxLength (unsigned int n, float *x, float *y, float *z, float *lx, float *ly, float *lz) {
   llog ("\n\tDetermining box length: ");
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
      llog ("\n\tShifting atoms by %f %f %f", off[X], off[Y], off[Z]);
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
   llog ("\n\tBox dimensions: %f %f %f %f %f %f", min[X], max[X], min[Y], max[Y], min[Z], max[Z]);
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
      llog ("%ld ", ptr->particle);
      ptr = ptr->next;
   }
}

/* --------------------------------------------------------------------------------------------- */

#define CoordinateToCell(x,y,z) (int) (    floor (x*CellLengthI[X]) +    floor (y*CellLengthI[Y])*NCellSide[X] +    floor (z*CellLengthI[Z])*NCellSide[X]*NCellSide[Y]   )


int *BuildLinkedList (int n, float *x, float *y, float *z, float RcNeighbor) {
   llog ("\nBuilding cell list...");
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

   llog ("\n\tNcells: %d %d %d -> %d    CellSize: %f   CellLength: %f %f %f", 
      NCellSide[X], NCellSide[Y], NCellSide[Z], NCell, CellSize, CellLength[X], CellLength[Y], CellLength[Z]);

   CellList = (t_linkedlist **) calloc (NCell, sizeof (t_linkedlist));
   llog ("\n\tMemory: %fMB", (sizeof (t_linkedlist)*NCell+n*sizeof (unsigned int))/1024./1024.);
   
   for (i = 0; i < n; i++) 
      LinkedListInsert (i, CellList + CoordinateToCell (x[i], y[i], z[i]) );

   return NULL;
}

/* --------------------------------------------------------------------------------------------- */

int *BuildNeighborList (int n, float *x, float *y, float *z, float RcNeighbor) {
   llog ("\nBuilding neighbor list...");
   return BuildLinkedList (n, x, y, z, RcNeighbor);
}

/* --------------------------------------------------------------------------------------------- */

void MoleculeListNeighborLinkedList (int Atoms, float *x, float *y, float *z, float RcClusterSq) {
   llog ("\nBegin clustering...");
   double start = clock();
   
   float z0=0.0;
   float RcCluster = sqrt (RcClusterSq);
   int i, j, k, Lk, Size, s, n;

   float xCell[3];
   int m[3];
   t_linkedlist *ptr;

   t_linkedlist **List = (t_linkedlist **) calloc (Atoms, sizeof (t_linkedlist));
   int *Spectrum = (int *) calloc (Atoms+1, sizeof (int));
   if (!List || !Spectrum)
      exit (EMEM);
   llog ("\n\t-> %fMB memory\n", (double) (sizeof (t_linkedlist)*Atoms + sizeof (int)*(Atoms+1))/1024./1024. );
   
   for (i = 0; i < Atoms; i++) {
      List[i]->particle = i;
      List[i]->next = List[i];
   }
   llog ("\n\tstart");
   #define NOT_YET_CLUSTERED(i) (List[i] == List[i]->next)
   for (i = 0; i < Atoms; i++) {
      PROGRESS(i);
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
                           if (RSQ(j,k) <= RcClusterSq) { 
                              List[k]->next = List[j]->next;
                              List[j]->next = List[k]->next;
                           }
                        }
                     }

                  }
               }
            }
            // set j
            j = List[j]->next->particle;

            PROGRESS(j);
         }  while (i != j);
      }  /* if i */
   }  /* for i */
   
   double stop = clock();
   llog ("\n\tfinished in %.1f seconds", (stop-start)/CLK);
/*
   for (i = 0; i < Atoms+1; i++) 
      Spectrum[i] = 0;
   for (i = 0; i < Atoms; i++) {
      if (List[i] >= 0) {
         Size = 1;  j = List[i];  List[i] = -1;
         while (j != i) { Size++;   k = List[j];  List[j] = -1;  j = k; }
         //Bestimmung der Groesse des Clusters im Startteilchen i
         Spectrum[Size] += 1;
      }  
   }  
   printf("\nCoherent Clusters below %4.4f A\n", RcClusterSq);
   printf("Size Times Atoms\n");
   for (i = 0; i < Atoms+1; i++) // Atoms + 1 because 0 is also a possible Count   
   {
      if ((s = Spectrum[i]) != 0) 
         printf ("%6d %6d %6d\n", i, s, i*s);
   }
*/
   free (List);
   free (Spectrum);
   #undef NOT_YET_CLUSTERED
}

/* --------------------------------------------------------------------------------------------- */

void MoleculeListNeighbor (int Atoms, float *x, float *y, float *z, float RcClusterSq) {
   llog ("\nBegin clustering...");
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
   llog ("\n\t-> %fMB memory\n", (double) (sizeof (INT)*(Atoms*2+1))/1024./1024. );

   for (i = 0; i < Atoms; i++)  
      List[i] = i;

   #define NOT_YET_CLUSTERED(i) (i ==  List[i])
   for (i = 0; i < Atoms; i++) {
      PROGRESS(i);
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

            PROGRESS(j);
         }  while (i != j);
      }  /* if i */
   }  /* for i */
   
   double stop = clock();
   llog ("\n\tfinished in %.1f seconds", (stop-start)/CLK);

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
   llog ("\nBegin clustering...");
   double start = clock();
   
   FLOAT     Xj, Yj, Zj, RSq, z0=0.0;

   INT        i, j, k, Lk, Size, *List,*Spectrum,s,n;

   List = (INT *) malloc (sizeof (INT) * Atoms);
   Spectrum = (INT *) malloc (sizeof (INT) * Atoms+1);
   if (!List || !Spectrum)
      exit (EMEM);
   llog ("\n\t-> %fMB memory\n", (double) (sizeof (INT)*(Atoms*2+1))/1024./1024. );

   for (i = 0; i < Atoms; i++)  
      List[i] = i;

   for (i = 0; i < Atoms-1; i++) {
      PROGRESS(i);
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
            PROGRESS(j);
         }  while (i != j);
      }  /* if i */
   }  /* for i */
   
   double stop = clock();
   llog ("\n\tfinished in %.1f seconds", (stop-start)/CLK);

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

void ParseArgs (int argc, char *argv[], float *rcut, unsigned int *Neighborlist, float *RcNeighbor, 
      unsigned int *LinkedList, unsigned int *FileFlag, char **filename) {
   int opt;
   while ((opt = getopt(argc, argv, "h:r:nlf:XYZ")) != -1) {
      switch (opt) {
         case 'h': printf (USAGE); break;
         case 'r': *rcut = (float) atof (optarg); break;
         case 'n': *Neighborlist = TRUE; /**RcNeighbor = (float) atof (optarg); */break;
         case 'l': *LinkedList = TRUE; break;
         case 'f': *filename = optarg; *FileFlag = TRUE; break;
         case 'X': Periodic[X] = TRUE; break;
         case 'Y': Periodic[Y] = TRUE; break;
         case 'Z': Periodic[Z] = TRUE; break;
         default:  printf (USAGE); break;
      }
   }
   llog ("\n\tPeriodicity: %d %d %d", Periodic[X], Periodic[Y], Periodic[Z]);
}

/* --------------------------------------------------------------------------------------------- */

void ReadStdin (unsigned int *n, float **x, float **y, float **z) {
   unsigned int i, ret;
   double start = clock();
   llog ("\nReading from stdin...");
   llog ("\n\tExpecting (x, y, z)");
   *n = 0;
   *x = NULL;
   *y = NULL;
   *z = NULL;
   float tmp;
   do {
      *n += 1;
      PROGRESS(*n);
      i = *n-1;
      *x = (float *) realloc (*x, *n*sizeof (float));
      *y = (float *) realloc (*y, *n*sizeof (float));
      *z = (float *) realloc (*z, *n*sizeof (float));
      if (!*x || !*y || !*z)
         exit (EMEM);
      ret = scanf ("%f %f %f\n", &((*x)[i]), &((*y)[i]), &((*z)[i]));
      fflush(stdin);
   } while (ret != EOF);
   *n -= 1;
   double stop = clock();
   llog ("\n\tread %d atoms -> memory %gMB in %.1f seconds", *n, (3* *n *sizeof (float))/(1024.*1024.), (stop-start)/CLK);
}

/* --------------------------------------------------------------------------------------------- */

void ReadFile (char *filename, unsigned int *n, float **x, float **y, float **z) {
   unsigned int i, ret;
   double start = clock ();
   llog ("\nReading from %s...", filename);
   llog ("\n\tExpecting (x, y, z)");
   *n = 0;
   *x = NULL;
   *y = NULL;
   *z = NULL;
   float tmp;
   FILE *in = fopen (filename, "r");

   char command[2048];
   sprintf (command, "cat %s | awk '{n++}END{print n}'", filename);
   llog ("\n\tUsing %s", command);
   FILE *pipe = popen (command, "r");
   if (!pipe) 
      exit (EFILE);
   char pos[1024];
   fgets (pos, sizeof(pos), pipe);
   fclose (pipe);
   unsigned int lines = (unsigned int) atoi (pos);
   llog ("\n\tWilling to read %d lines", lines);

   if (!in)
      exit (EFILE);
   do {
      *n += 1;
      PROGRESS(*n);
      i = *n-1;
      *x = (float *) realloc (*x, *n*sizeof (float));
      *y = (float *) realloc (*y, *n*sizeof (float));
      *z = (float *) realloc (*z, *n*sizeof (float));
      if (!*x || !*y || !*z)
         exit (EMEM);
//      ret = fscanf (in, "%f %f %f %f\n", &((*x)[i]), &((*y)[i]), &((*z)[i]), &tmp);
      ret = fscanf (in, "%f %f %f\n", &((*x)[i]), &((*y)[i]), &((*z)[i]));
   } while (*n-1 < lines);//ret != EOF);
   llog ("\nbreak with %d\n", ret);
   *n -= 1;
   fclose (in);
   double stop = clock();
   llog ("\n\tread %d atoms -> memory %gMB in %.1f seconds", *n, (3* *n *sizeof (float))/(1024.*1024.), (stop-start)/CLK);
}

/* --------------------------------------------------------------------------------------------- */

void ReadHDF (char *filename, unsigned int *n, float **x, float **y, float **z) {
   llog ("\nReading hdf file %s", filename);
   char *src = "/Particles/All";
   hid_t hdf = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
   hid_t particles = H5Dopen (hdf, src);
   hid_t type = H5Dget_type (particles);
   hid_t size = H5Dget_space (particles);
   hid_t plist = H5Dget_create_plist (particles);
   
   hsize_t nfields, nrecords;
   H5TBget_table_info ( hdf, src, &nfields, &nrecords );
   llog ("\thas %d fields with %d records\n", nfields, nrecords);

   char **field_names = (char **)calloc (nfields, sizeof (char*));
   int i;
   for (i = 0; i < nfields; i++) {
         field_names[i] = (char *)calloc (nfields, sizeof (char[256]));
   }
   size_t *field_sizes = (size_t *)calloc (nfields, sizeof (size_t)),
          *field_offsets = (size_t *)calloc (nfields, sizeof (size_t)),
          *type_size = (size_t *)calloc (nfields, sizeof (size_t));

   H5TBget_field_info (hdf, src, field_names, field_sizes, field_offsets, type_size);
   for (i = 0; i < nfields; i++)
      llog ("\t-> %s {%d}\n", field_names[i], field_sizes[i], field_offsets[i], type_size[i]);
  
   size_t total = 0;
   for (i = 0; i < nfields; i++)
      total += field_sizes[i];
   llog ("\ntotal size:%d\n", total);

   void *dst_buf = NULL;//calloc (10, total);
   
   hsize_t start = 0, nn = 1;
   char *field = "x";
   int fidx = -1;
   for (i = 0; i < nfields; i++) {
      if (strcmp (field_names[i], field) == 0) {
         fidx = i;
      }
   }
//   H5TBread_records (hdf, src, start, nn, type_size[0], field_offsets, field_sizes, dst_buf);
   dst_buf = calloc (nn, sizeof(float)); 
   H5TBread_fields_name (hdf, src, field, start, 0, type_size[0],  field_offsets, field_sizes, dst_buf);

   void *ptr = dst_buf;
//   for (i = 0; i < fidx; i++)
  //    ptr += field_sizes[i];

   float *variable1 = (float *) ptr;

   llog ("\ndone\n");
   //printf ("\nok data:\n");
   //for (i = 0; i < nfields; i++)
   //   printf ("%f ", dst_buf[i]);
   hid_t x1 = H5Tcreate (H5T_COMPOUND, sizeof (float));
//   H5Tinsert (x1, "x", HOFFSET(x1,
//   H5Dread (particles,  
   
//   H5TBread_table (hdf, "/Particles/All", type_size[0], field_offsets, field_sizes, dst_buf);

//   HD5read (particles, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dst_buf);
   
   H5Dclose (particles);
   H5Fclose (hdf);
   exit (0);
}

/* --------------------------------------------------------------------------------------------- */

int main (int argc, char *argv[]) {
   unsigned int n;
   float *x, *y, *z;

   float RcClusterSq, RcCluster = 0.;
   float RcNeighbor = 0;
   unsigned int Neighborlist = FALSE,
      LinkedList = FALSE;

   unsigned int FileFlag = FALSE;
   char *filename[1024];
   Periodic[X] = Periodic[Y] = Periodic [Z] = FALSE;

   ParseArgs (argc, argv, &RcCluster, &Neighborlist, &RcNeighbor, &LinkedList, &FileFlag, filename);

   if (RcCluster == 0.) 
      exit (ESYN);
   else
      RcClusterSq = SQ(RcCluster);
   
   if (!FileFlag)
      ReadStdin (&n, &x, &y, &z);
   else
      ReadHDF (*filename, &n, &x, &y, &z);
   
   DetermineBoxLength (n, x, y, z, &(Length[X]), &(Length[Y]), &(Length[Z]));

   if (Neighborlist) {
      RcNeighbor = RcCluster;
      BuildNeighborList (n, x, y, z, RcNeighbor);
      if (LinkedList) 
         MoleculeListNeighborLinkedList (n, x, y, z, RcClusterSq);
      else 
         MoleculeListNeighbor (n, x, y, z, RcClusterSq);
   } else 
      MoleculeList (n, x, y, z, RcClusterSq);

   return EOK;
}
