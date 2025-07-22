#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <stdbool.h>
#define RAVE_N_ATOMS 500
#define LARGE_ENERGY 1E20

/******************** Prototypes *********************************/
struct Area {
              double sas, ctc, pc_sas, pc_ctc, sas_self, ctc_self, sas_res, ctc_res;
            } Area;

struct System_area {
               struct Area system;
               struct Area *atoms;
               struct Area *residues;
               struct Area *chains;
               struct Area *segments;
               struct Area *models;
             } System_area;

struct At_srf
{
int n_isrf; // numero di punti di superficie
int *isrf;
};

struct Srf_pt
{
double *r;
double *vec;
double a;
int at, atn;
int type;
};

struct Srf
{
int alloc_srf_pt;
int n_srf_pt_tot;
struct Srf_pt *srf_pt;
struct At_srf *at_srf;
};

struct System { 
               int n_atoms;
               struct Atom *atoms;
               int n_residues;
               struct Residue *residues;
               int n_chains;
               struct Chain *chains;
               int n_segments;
               struct Segment *segments;
               int n_models;
               struct Model *models;
              }  System;

struct Atom_grid {
                  int **atom_node;
                  int PBC_ON; /*periodic boundary conditions 1 or 0 */
                  double pbcx, pbcy, pbcz;
                  int *n_atom_node;
                  int grid_X, grid_Y, grid_Z, grid_size, max_atom_node;
                  double x_min, y_min, z_min, x_max, y_max, z_max;
                  double mesh;
                  double dx,dy,dz;
                 } Atom_grid;

struct Neighbour {
                  int max_n;             /* max n neighbours */
                  int n_neighbours;             /* number of neighbours */
                  int *list; 
                  double *d;
                 } Neighbour; /* list of atom index of neighbours */


struct Atom {
/************** fields in PDB ********************/ 
		char at_name[5]; 
		char alt_loc;
		char res_name[4];
		char chain; 
                char element[3]; 
		int model;
                int at_n;
		int res_n;
		char res_ins; 
                double coor[3]; /* queste sostituiranno x, y, z */
                double ref[3];
		double occ,temp;
		char segid[5];
     		char pdb_chrg[3]; /* stringa con la carica formale in PDB */
/************* additional features **********************/
    		double mass;
    		char chem_group[5];  /* chemical group e.g. peptide methyl.... */
    		char group[5];  /* if sidechain or backbone or other */
    		double radius;
    		double charge;
    		int atom_type;
		int i;
	    } Atom;

struct Trj {
           int nf, naxf;  // number of frames, number of atoms per frame
           double **coor; 
           } Trj;

struct Residue { char res_name[4];
                 char chain;
                 char res_ins;
                 int res_n;
                 int model;
                 char seg_name[5];
                 int beg, end ; /*begin and end atoms of sorted atom list*/
                 int n_alt_loc;  
                 int prev, next;
                 double CA[3], CG[3], CM[3]; 
                 int res_type;
               } Residue;

struct Chain { char ch_name;
               int beg, end;
               int model;
               char seg_name[5];
             } Chain;


struct Segment { char seg_name[5];
               int beg, end;
               int model;
               }  Segment;

struct Model { int model_n;
               int beg, end;
               }  Model;

FILE *file_open(char *fname, char *acc);

void make_system(struct System *system, struct Trj *trj);

/* make the different between two double vectors r1 and r2. The result is in v*/
void diffv(double *v, double *r2, double *r1);
/* make the sum of two double vectors r1 and r2. The result is in v*/
void sumv(double *v, double *r2, double *r1);
/*multiply a vector by a scalar*/
void cmul(double *v, double a, double *v0);
void scale(double *v, double a);
/* funzione per il calcolo della distanza tra due vettori r1 e r2*/
double distv(double *r1, double *r2);
double distv_pbc(double *r1, double *r2, double *s);
double dist(double x1, double y1, double z1, double x2, double y2, double z2);

/* calcolo dell'angolo tra tre vettori*/
double anglev(double *x1, double *x2, double *x3);
double angle(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3);

/* calcolo dell'angolo (pseudo??)torsionale tra i vettori x1, x2, x3, x4*/
double torsionv(double *x1, double *x2, double *x3, double *x4);
double torsion(double x1, double y1, double z1, 
																				double x2, double y2, double z2,
                    double x3, double y3, double z3, 
	                   double x4, double y4, double z4);
	                    
/* funzione per il calcolo del prodotto scalare di due vettori x1 e x2*/
double dotv(double *x1, double *x2);
double dot(double x1, double y1, double z1, double x2, double y2, double z2);
double dotvn(double *x1, double *x2,int n);
double modv(double *x1);
double modvn(double *x1, int n);

/* calcola il prodotto vettoriale di due vettori reali x1 e x2, in x*/
void vecv(double *x, double *x1, double *x2);
void vec(double *x, double *y, double *z,
         double x1, double y1, double z1, double x2, double y2, double z2);

void make_system(struct System *system, struct Trj *trj);


/******************* includes  *********************/
void hpsort(struct Atom *ra, int n);
void hpsort_trj(struct Atom *ra, int n, struct Trj *trj);

/* Prototypes */
void read_atom_pdb(char *buf, struct Atom *atom);

void read_PQR_atoms(FILE *fp1, int *n_atoms, struct Atom **atoms, struct Trj *trj);

void read_atom_pqr(char *buf, struct Atom *atom);

void write_PDB_atoms(FILE *fp, int n_atoms, struct Atom *atoms, struct Trj trj);

void write_PQR_atoms(FILE *fp1, int n_atoms, struct Atom *atoms, struct Trj trj);

void write_atom_pdb(char *buf, struct Atom atom);

void write_atom_pqr(char *buf, struct Atom atom);

void copy_atom(struct Atom *atom_to, struct Atom atom_from);

void get_grid_parameters(struct System sistema, double probe_radius,
                double cutoff, struct Atom_grid *atom_grid, char type);

void print_grid_info(struct Atom_grid atom_grid);

void atoms_on_grid(struct Atom_grid *atom_grid, struct System sistema);

void grid_neighbour(struct Atom_grid atom_grid, 
                struct System sistema, double cutoff, struct Neighbour *neighbour, char type);

void system_area_alloc(struct System system, struct System_area *system_area);

void system_area_free(struct System_area *system_area);

void sort_neighbours(struct Neighbour *neighbour, int beg, int end);

void sort_neighbours_by_dist(struct Neighbour *neighbour, int beg, int end);

struct Neighbour *neighbour_alloc(int n, int max_nei);

void neighbour_realloc(struct Neighbour *neighbours, int n);

void neighbour_free(struct Neighbour *neighbours, int n);

void atom_grid_alloc(struct Atom_grid *atom_grid);

void atom_grid_realloc(struct Atom_grid *atom_grid);

void atom_grid_free(struct Atom_grid *atom_grid);

void srf_free(struct Srf *srf, int n_atoms);

int is_aa(struct Residue residue);

int is_h_ok(struct Residue residue);

void eliminate_alt_loc(struct System *sorted_system, struct Trj *trj);

void guess_mass(struct System *system);

int next_res(struct System system, int m, int n);

int cmp_atoms(const void *p1, const void *p2);


/* confronta due atomi*/
int cmp_atoms(const void *p1, const void *p2)
{
	struct Atom A_atom, B_atom;
	int check = 0 ; 

	A_atom = *((struct Atom *)p1); 
	B_atom = *((struct Atom *)p2); 


/**** confronto    if (numbers[min] <= numbers[mid]) *****/

	if( A_atom.model < B_atom.model) check = -1; 
	else if (A_atom.model == B_atom.model) 
    {
        check = 0;
    	if( strcmp(A_atom.segid,B_atom.segid) < 0) check = -1; 
    	else if (!strcmp(A_atom.segid,B_atom.segid)) 
        {
        	check = 0;
           	if( (int) A_atom.chain < (int) B_atom.chain ) check = -1;
           	else if( A_atom.chain ==  B_atom.chain )
            {
            	check = 0;
                if( A_atom.res_n <  B_atom.res_n ) check = -1;
                else if(A_atom.res_n ==  B_atom.res_n )
                {
                	check = 0;
                    if((int) A_atom.res_ins < (int) B_atom.res_ins ) check = -1;
                    else if( A_atom.res_ins ==  B_atom.res_ins )
                    {
                    	check = 0;
                        if( strcmp(A_atom.at_name,B_atom.at_name) < 0) check = -1; 
                        else if (!strcmp(A_atom.at_name,B_atom.at_name)) 
                        {
                        check = 0;
                        if( (int) A_atom.alt_loc < (int) B_atom.alt_loc ) check = -1;
                        else if( A_atom.alt_loc ==  B_atom.alt_loc ) check = 0;
                        else check = 1;
			}
                        else check = 1;
                    }
                        else check = 1;
                }
                        else check = 1;
            }
                        else check = 1;
        }
                        else check = 1;
	}
        else check = 1;
	return check;

}

void eliminate_alt_loc(struct System *system, struct Trj *trj)
{
int i,j,k,l;
char buf[256];
j = 1;
for (i=1; i< (*system).n_atoms; i++)   
{
/*
write_atom_pdb(buf,(*system).atoms[i]);
printf("%s\n",buf);
*/
if(
strcmp((*system).atoms[i-1].at_name,(*system).atoms[i].at_name) ||
strcmp((*system).atoms[i-1].res_name,(*system).atoms[i].res_name) ||
strcmp((*system).atoms[i-1].segid,(*system).atoms[i].segid) ||
((*system).atoms[i-1].chain != (*system).atoms[i].chain) ||
((*system).atoms[i-1].model != (*system).atoms[i].model) ||
((*system).atoms[i-1].res_n != (*system).atoms[i].res_n) ||
((*system).atoms[i-1].res_ins != (*system).atoms[i].res_ins)
)
{
//printf("sono qui %i/%i\n",i,(*system).n_atoms);
if(i != j)
{
copy_atom(&((*system).atoms[j]), ((*system).atoms[i]));
for(k = 0; k < (*trj).nf; k++)
for(l = 0; l<3; l++)
(*trj).coor[k * (*trj).naxf + j][l] = (*trj).coor[k * (*trj).naxf + i][l];
}
/* toglie la alt_loc */
(*system).atoms[j].alt_loc = ' ';
j++;
//if((i+1) == (*system).n_atoms-1) 
//copy_atom(&((*system).atoms[j]), ((*system).atoms[i+1]));
}
else 
{
// Commentate per evitare eccessivo output
//printf("Elimino l'atomo:\n");
//write_atom_pdb(buf,(*system).atoms[i]);
//printf("%s\n",buf);
//printf("perche' c'e' gia':\n");
//write_atom_pdb(buf,(*system).atoms[i-1]);
//printf("%s\n",buf);
}

}
if ((*system).n_atoms > 0)
{
(*system).n_atoms = j;
for(k = 0; k < (*trj).nf; k++)
for(j = 0; j< (*system).n_atoms; j++)
for(l = 0; l<3; l++)
(*trj).coor[k * (*system).n_atoms + j][l] = (*trj).coor[k * (*trj).naxf + j][l];
(*trj).naxf = (*system).n_atoms;
}
}


/*funzione per il controllo del nome di un residuo*/ 
int is_aa(struct Residue residue)
{
int p;
if(!strncmp(residue.res_name, "ALA", 3) ||
   !strncmp(residue.res_name, "CYS", 3) || 
   !strncmp(residue.res_name, "CYM", 3) || 
   !strncmp(residue.res_name, "CYX", 3) || 
   !strncmp(residue.res_name, "CSS", 3) || 
   !strncmp(residue.res_name, "ASP", 3) || 
   !strncmp(residue.res_name, "GLU", 3) || 
   !strncmp(residue.res_name, "PHE", 3) || 
   !strncmp(residue.res_name, "GLY", 3) || 
   !strncmp(residue.res_name, "HIS", 3) || 
   !strncmp(residue.res_name, "HID", 3) || 
   !strncmp(residue.res_name, "HIE", 3) || 
   !strncmp(residue.res_name, "HIP", 3) || 
   !strncmp(residue.res_name, "HSD", 3) || 
   !strncmp(residue.res_name, "HSE", 3) || 
   !strncmp(residue.res_name, "HSP", 3) || 
   !strncmp(residue.res_name, "ILE", 3) || 
   !strncmp(residue.res_name, "LYS", 3) || 
   !strncmp(residue.res_name, "LEU", 3) || 
   !strncmp(residue.res_name, "MET", 3) || 
   !strncmp(residue.res_name, "ASN", 3) || 
   !strncmp(residue.res_name, "PRO", 3) || 
   !strncmp(residue.res_name, "GLN", 3) || 
   !strncmp(residue.res_name, "ARG", 3) || 
   !strncmp(residue.res_name, "SER", 3) || 
   !strncmp(residue.res_name, "THR", 3) || 
   !strncmp(residue.res_name, "VAL", 3) || 
   !strncmp(residue.res_name, "TRP", 3) || 
   !strncmp(residue.res_name, "TYR", 3) ) p = 1;
else p = 0;
return p;
}

int next_res(struct System system, int m, int n)
{
	int k,found1, found2, p;
	struct Residue *residue1 = &(system.residues[m]);
	struct Residue *residue2 = &(system.residues[n]);
        double *v1, *v2;
	found1 = found2 = 0;
	if(is_aa(*residue1) && is_aa(*residue2) )
	{
		k=(*residue1).beg;
		while((k <= (*residue1).end) && (!found1)) {
			if(!strcmp(system.atoms[k].at_name,"CA"))
                        {       v1 = system.atoms[k].coor;
				found1 = 1;
                        }
			k++;
		}
		k=(*residue2).beg;
		while((k <= (*residue2).end) && (!found2)) {
			if(!strcmp(system.atoms[k].at_name,"CA"))
                        {       v2 = system.atoms[k].coor;
				found2 = 1;
                        }
			k++;
		}
               if( found1 && found2 && (distv(v1,v2) < 4.0)) p = 1; 
               else p = 0;
	}
	return p;
}

/* creazione ed inizializzazione strutture di sistema*/
void make_system(struct System *system, struct Trj *trj)
{
	int n_atoms = system->n_atoms;
	struct Atom *atoms = system->atoms;
	struct Residue *residues = system->residues; 
	int *p_n_residues = &(system->n_residues); 
	int *p_n_chains = &(system->n_chains); 
	struct Chain *chains = system->chains; 
	int *p_n_segments  = &(system->n_segments); 
	struct Segment *segments = system->segments;
	int *p_n_models = &(system->n_models); 
	struct Model *models = system->models; 

	int imodel=-1, isegment=-1, ichain=-1, iresidue=-1;
	char res_ins = '*';
	char chain = '*';
	char segid[5] = "*****";
	int model = -1;

	int i;
	int res_n;

//        if( system->n_atoms < 100000)
//        {
//        printf("sorting %i atoms using quicksort\n", system->n_atoms);  
//	qsort(system->atoms, system->n_atoms, sizeof(struct Atom), &cmp_atoms);
//        }
//        else
//      {

        printf("sorting %i atoms using heapsort\n", system->n_atoms);  
        for(i=0;i<system->n_atoms;i++) system->atoms[i].i = i;
	hpsort_trj(system->atoms, system->n_atoms,trj);
//	hpsort(system->atoms, system->n_atoms);
//        }

        eliminate_alt_loc(system, trj);
        printf("after eliminate_alt_loc %i atoms left\n", system->n_atoms);  
	n_atoms = system->n_atoms;
// qui conta models, segments, chains, residue,
        model = -1;
	for (i=0; i<n_atoms; i++)
		if(atoms[i].model != model)
		{
//                        printf("atom %i: model: %i vs. %i\n", i, atoms[i].model, model);
			imodel++;
			isegment++;
			ichain++;
			iresidue++;
                        model=atoms[i].model;
                        strcpy(segid,atoms[i].segid);
                        chain=atoms[i].chain;
                        res_n=atoms[i].res_n;
                        res_ins=atoms[i].res_ins;
                }
		else if(strcmp(atoms[i].segid,segid))
			{
 //                       printf("atom %i: segid: %s vs. %s\n", i, atoms[i].segid, segid);
			isegment++;
			ichain++;
			iresidue++;
                        strcpy(segid,atoms[i].segid);
                        chain=atoms[i].chain;
                        res_n=atoms[i].res_n;
                        res_ins=atoms[i].res_ins;
                        }
		else if(atoms[i].chain != chain) 
                        {
			ichain++;
			iresidue++;
                        chain=atoms[i].chain;
                        res_n=atoms[i].res_n;
                        res_ins=atoms[i].res_ins;
                        }
		else if(atoms[i].res_n != res_n) 
                       {
			iresidue++;
                        res_n=atoms[i].res_n;
                        res_ins=atoms[i].res_ins;
                       }
		else if(atoms[i].res_ins != res_ins) 
                       {
			iresidue++;
                        res_ins=atoms[i].res_ins;
                       }
    imodel++;
    isegment++;
    ichain++;
    iresidue++;
    system->models = (struct Model *) calloc((size_t) imodel, sizeof(Model));
    if(system->models==NULL) {printf("Could not allocate memory for System.models... Exiting...\n"); exit(0);}
    system->segments = (struct Segment *) calloc((size_t) isegment, sizeof(Segment));
    if(system->segments==NULL) {printf("Could not allocate memory for System.segments... Exiting...\n"); exit(0);}
    system->chains = (struct Chain *)  calloc((size_t) ichain, sizeof(Chain));
    if(system->chains==NULL) {printf("Could not allocate memory for System.chains... Exiting...\n"); exit(0);}
    system->residues = (struct Residue *) calloc((size_t) iresidue, sizeof(Residue));
    if(system->residues==NULL) {printf("Could not allocate memory for System.residues... Exiting...\n"); exit(0);}

//        printf("%i %i %i %i\n", imodel, isegment, ichain, iresidue);
	imodel=-1, isegment=-1, ichain=-1, iresidue=-1;
        model = -1;
	for (i=0; i<n_atoms; i++)
	{
		if(atoms[i].model != model)
		{
			imodel++;
			isegment++;
			ichain++;
			iresidue++;
//        printf("%i %i %i %i\n", imodel, isegment, ichain, iresidue);
 
			strcpy(system->residues[iresidue].res_name, atoms[i].res_name);
			system->residues[iresidue].res_n = atoms[i].res_n;
			system->residues[iresidue].res_ins = atoms[i].res_ins;
			system->residues[iresidue].chain = atoms[i].chain;
			system->residues[iresidue].model = atoms[i].model;
			strcpy(system->residues[iresidue].seg_name, atoms[i].segid);
			system->residues[iresidue].beg = i;
			if (iresidue != 0) system->residues[iresidue - 1].end = i-1;
        	chain = atoms[i].chain; 
        	res_n = atoms[i].res_n; 
        	res_ins = atoms[i].res_ins; 
  
/*		if (sscanf("%s", atoms[i].segid) == 0) 
            	strcpy(segid,"     "); 
        	else */
            	strcpy(segid, atoms[i].segid); 
			model=atoms[i].model;        

			system->chains[ichain].ch_name = chain;
			system->chains[ichain].beg = i;
                 	system->chains[ichain].model=atoms[i].model;
			strcpy(system->chains[ichain].seg_name, atoms[i].segid);
			if (ichain != 0) system->chains[ichain - 1].end = i-1;

			strcpy(system->segments[isegment].seg_name, segid);
        	system->segments[isegment].model=atoms[i].model;
        	system->segments[isegment].beg = i;
 			if (isegment != 0) system->segments[isegment - 1].end = i-1;

        	system->models[imodel].beg=i;
        	system->models[imodel].model_n=model;
        	if (imodel != 0) system->models[imodel - 1].end =i-1;
        //printf("Sono qui: model: %i prev: %i\n", atoms[i].model,model);
        //printf("1: %i %s\n", i, (*system).atoms[i].at_name);
		}
		else if(strcmp(atoms[i].segid,segid))
			{
				isegment++;
  				ichain++;
  				iresidue++;

				strcpy(system->residues[iresidue].res_name, atoms[i].res_name);
				system->residues[iresidue].res_n = atoms[i].res_n;
				system->residues[iresidue].res_ins = atoms[i].res_ins;
				system->residues[iresidue].chain = atoms[i].chain;
				system->residues[iresidue].model = atoms[i].model;
				strcpy(system->residues[iresidue].seg_name, atoms[i].segid);
				system->residues[iresidue].beg = i;
			/*	residues[iresidue].res_type = restyp(atoms[i].res_name); */
				if (iresidue != 0) system->residues[iresidue - 1].end = i-1;
			    chain = atoms[i].chain; 
    		    res_n = atoms[i].res_n; 
       		 	res_ins = atoms[i].res_ins; 
/*       	 		if (sscanf("%s", atoms[i].segid) == 0) 
           		 	strcpy(segid,"     "); 
        		else */
            		strcpy(segid, atoms[i].segid); 
        
				system->chains[ichain].ch_name = chain;
				system->chains[ichain].beg = i;
        		system->chains[ichain].model=atoms[i].model;
				strcpy(system->chains[ichain].seg_name, atoms[i].segid);
				if (ichain != 0) system->chains[ichain - 1].end = i-1;

	      		strcpy(system->segments[isegment].seg_name, segid);
    	    	system->segments[isegment].model=atoms[i].model;
       	 		system->segments[isegment].beg = i;
 				if (isegment != 0) system->segments[isegment - 1].end = i-1;
        //printf("Sono qui: segid: %s prev: %s\n", atoms[i].segid,segid);
        //printf("2: %i %s\n", i, (*system).atoms[i].at_name);

			}
			else if(atoms[i].chain != chain) 
				{
					ichain++;
					iresidue++;

					strcpy(system->residues[iresidue].res_name, atoms[i].res_name);
					system->residues[iresidue].res_n = atoms[i].res_n;
					system->residues[iresidue].res_ins = atoms[i].res_ins;
					system->residues[iresidue].chain = atoms[i].chain;
					system->residues[iresidue].model = atoms[i].model;
					strcpy(system->residues[iresidue].seg_name, atoms[i].segid);
					system->residues[iresidue].beg = i;
				/*	residues[iresidue].res_type = restyp(atoms[i].res_name); */
					if (iresidue != 0) system->residues[iresidue - 1].end = i-1;
        			chain = atoms[i].chain; 
       		 		res_n = atoms[i].res_n; 
        			res_ins = atoms[i].res_ins; 
					system->chains[ichain].ch_name = chain;
					system->chains[ichain].beg = i;
        			system->chains[ichain].model=atoms[i].model;
					strcpy(system->chains[ichain].seg_name, atoms[i].segid);
					if (ichain != 0) system->chains[ichain - 1].end = i-1;
        //printf("Sono qui: chain: %c prev: %c\n", atoms[i].chain,chain);
        //printf("3: %i %s\n", i, (*system).atoms[i].at_name);
				}
				else if(atoms[i].res_n != res_n) 
					{
					iresidue++;

					strcpy(system->residues[iresidue].res_name, atoms[i].res_name);
					res_n = atoms[i].res_n;
					res_ins = atoms[i].res_ins;
					system->residues[iresidue].res_n = res_n;
					system->residues[iresidue].res_ins = res_ins;
					system->residues[iresidue].chain = chain;
					system->residues[iresidue].model = atoms[i].model;
					strcpy(system->residues[iresidue].seg_name, atoms[i].segid);
					system->residues[iresidue].beg = i;
				/*	residues[iresidue].res_type = restyp(atoms[i].res_name); */
					if (iresidue != 0) system->residues[iresidue - 1].end = i-1;
        //printf("Sono qui: res_n: %i prev: %i\n", atoms[i].res_n,res_n);
        //printf("5: %i %s\n", i, (*system).atoms[i].at_name);
					}
					else if(atoms[i].res_ins != res_ins) 
					{
					iresidue++;
					strcpy(system->residues[iresidue].res_name, atoms[i].res_name);
					res_n = atoms[i].res_n;
					res_ins = atoms[i].res_ins;
					system->residues[iresidue].res_n = res_n;
					system->residues[iresidue].res_ins = res_ins;
					system->residues[iresidue].chain = chain;
					system->residues[iresidue].model = atoms[i].model;
					strcpy(system->residues[iresidue].seg_name, atoms[i].segid);
					system->residues[iresidue].beg = i;
					if (iresidue != 0) system->residues[iresidue - 1].end = i-1;
        //printf("Sono qui: res_ins: %c prev: %c\n", atoms[i].res_ins,res_ins);
        //printf("6: %i %s\n", i, (*system).atoms[i].at_name);
					}
	}
if(n_atoms != 0)
    {
    system->segments[isegment].end = i-1;
	system->chains[ichain].end = i-1;
	system->residues[iresidue].end = i-1;
    system->models[imodel].end = i-1;
    }
	*p_n_chains = ichain+1;
	*p_n_segments = isegment+1;
	*p_n_residues = iresidue+1;
	*p_n_models = imodel+1;

	printf("########################################################\n"); 
	printf("# SYSTEM:                                              #\n"); 
	printf("########################################################\n\n"); 
	printf("atoms = %8i, residues = %8i, chains = %8i, segments = %8i, models = %8i\n\n", n_atoms, *p_n_residues, *p_n_chains,*p_n_segments,*p_n_models  );
	printf("########################################################\n\n");
}
//void hpsort(struct Atom *ra, int n, int *cmp_atoms)
void hpsort(struct Atom *ra, int n)
{  
    int N, i, parent, child;  
    struct Atom rra;  
    N = n;
    i = n/2;
    for (;;) { 
        if (i > 0) { 
            i--;           
            copy_atom(&(rra),(ra[i]));
//             t = arr[i];    
        } else {     
            n--;           
            if (n == 0) return; 
            copy_atom(&(rra),(ra[n]));
//            t = arr[n];    
            copy_atom(&(ra[n]),(ra[0]));
//            arr[n] = arr[0]; 
        }  
  
        parent = i; 
        child = i*2 + 1; 
  
        while (child < n) {  
//            if (child + 1 < n  &&  arr[child + 1] > arr[child]) {  
            if (child + 1 < n  && (cmp_atoms(&(ra[child + 1]),&(ra[child])) > 0)) {
                child++; 
            }  
            if (cmp_atoms(&(ra[child]),&(rra)) > 0) {   
                copy_atom(&(ra[parent]),(ra[child]));
//                arr[parent] = arr[child]; 
                parent = child;  
                //child = parent*2-1; 
                child = parent*2+1;
            } else {  
                break; 
            }  
        }  
       copy_atom(&(ra[parent]),(rra));
//   arr[parent] = t;   
    }  
}
void hpsort_trj(struct Atom *ra, int n, struct Trj *trj) // qui da fare ancora la parte per la traiettoria
{  
    int N, i, parent, child, j, k;  
    struct Atom rra;  
    double **tmptrj;

    tmptrj = (double **) calloc((size_t) n,sizeof(double *));
    for(i = 0; i < n; i++) 
      tmptrj[i] = (double *) calloc((size_t) 3, sizeof(double));

//printf("sono qui 2\n");
// per la traiettoria do ad ogni atomo un indice, ordino prima gli atomi
// e poi con gli indici ordino la traiettoria
//    for(i = 0; i < n; i++) ra[i].i = i;
   N = n;
    i = n/2;
    for (;;) { 
        if (i > 0) { 
            i--;           
            copy_atom(&(rra),(ra[i]));
//             t = arr[i];    
        } else {     
            n--;           
            if (n == 0) goto final_ops;
            copy_atom(&(rra),(ra[n]));
//            t = arr[n];    
            copy_atom(&(ra[n]),(ra[0]));
//            arr[n] = arr[0]; 
        }  
  
        parent = i; 
        child = i*2 + 1; 
  
        while (child < n) {  
//            if (child + 1 < n  &&  arr[child + 1] > arr[child]) {  
            if (child + 1 < n  && (cmp_atoms(&(ra[child + 1]),&(ra[child])) > 0)) {
                child++; 
            }  
            if (cmp_atoms(&(ra[child]),&(rra)) > 0) {   
                copy_atom(&(ra[parent]),(ra[child]));
//                arr[parent] = arr[child]; 
                parent = child;  
                //child = parent*2-1; 
                child = parent*2+1;
            } else {  
                break; 
            }  
        }  
       copy_atom(&(ra[parent]),(rra));
//   arr[parent] = t;   
    }  
    //qui la traiettoria
final_ops:
       n = (*trj).naxf;
//       for(j = 0; j < n; j++)
//        printf("%i %i\n", j, ra[j].i);
    for(i = 0; i < (*trj).nf; i++)   
       {
       for(j = 0; j < n; j++)
         for(k = 0; k < 3; k++)
            tmptrj[j][k] = (*trj).coor[i * (*trj).naxf + ra[j].i][k];
       for(j = 0; j < n; j++)
         for(k = 0; k < 3; k++)
            (*trj).coor[i * (*trj).naxf + j][k] =  tmptrj[j][k]; 
       }
}
void cmul(double *v, double a, double *v0)
{
	v[0] = v0[0] * a;
	v[1] = v0[1] * a;
	v[2] = v0[2] * a;
}

void scale(double *v, double a)
{
	v[0] = v[0] * a;
	v[1] = v[1] * a;
	v[2] = v[2] * a;
}

void diffv(double *v, double *r2, double *r1)
{
	v[0] = r2[0] - r1[0];
	v[1] = r2[1] - r1[1];
	v[2] = r2[2] - r1[2];
}

/* v = r2 + r1*/
void sumv(double *v, double *r2, double *r1)
{
	v[0] = r2[0] + r1[0];
	v[1] = r2[1] + r1[1];
	v[2] = r2[2] + r1[2];
}

/* restituisce il prodotto scalare di due vettori ??*/
double dotv(double *x1, double *x2)
{       
    double d;
    d = x1[0]*x2[0] +  x1[1]*x2[1] + x1[2]*x2[2];
    return d;
}

double modv(double *x1)
{       
    double d;
    d = sqrt(x1[0]*x1[0] +  x1[1]*x1[1] + x1[2]*x1[2]);
    return d;
}

double dotvn(double *x1, double *x2, int n)
{       
    int i;
    double d=0.0;
    for(i=0; i<n; i++)
    d = d + x1[i]*x2[i];
    return d;
}

double modvn(double *x1, int n)
{       
    int i;
    double d=0.0;
    for(i=0; i<n; i++)
    d = d + x1[i]*x1[i];
    d = sqrt(d);
    return d;
}

/* prodotto vettoriale di x1 e x2*/
void vecv(double *x, double *x1, double *x2)
{
    x[0] = x1[1]*x2[2]-x1[2]*x2[1];
    x[1] = x1[2]*x2[0]-x1[0]*x2[2];
    x[2] = x1[0]*x2[1]-x1[1]*x2[0];
}
 
/* distanza (??) tra due vettori*/
double distv(double *r1, double *r2)
{
	double d,v[3];
	diffv(v,r2,r1);
	d = dotv(v, v );
	d =  sqrt(d);
	return d;
}

/* distanza (??) tra due vettori*/
double distv_pbc(double *r1, double *r2, double *s)
{
        int i; 
	double d,v[3];
	diffv(v,r2,r1);
        for(i = 0; i<3; i++)
          if(v[i] > s[i] / 2) v[i] = v[i] - s[i];
          else
          if(v[i] < -s[i] / 2) v[i] = v[i] + s[i];

	d = dotv(v, v );
	d =  sqrt(d);
	return d;
}

/* calcolo dell'angolo tra tre punti*/
double anglev(double *x1, double *x2, double *x3)
{
      
    double d21,d32,x21[3],x32[3], f, angle;
    diffv(x21, x2, x1);
    diffv(x32, x3, x2);

    d21 = sqrt(dotv(x21,x21));
    d32 = sqrt(dotv(x32,x32));
    f = -dotv(x21,x32)/(d21 * d32);

    angle = 180.0 * acos(f) / M_PI;
    return angle;
}

/* calcolo dell'angolo torsionale fra 4 punti*/
double torsionv(double *x1, double *x2, double *x3, double *x4)
{
        double d21,d32,d43,x21[3],x32[3],x43[3],bxc[3];
        double  ac,ab,bc,abxc, t;
        int i;


        diffv(x21, x2, x1);
        diffv(x32, x3, x2);
        diffv(x43, x4, x3);
        d21 = sqrt(dotv(x21,x21));
        d32 = sqrt(dotv(x32,x32));
        d43 = sqrt(dotv(x43,x43));

        for(i = 0; i < 3; i++)
        {
         x21[i] = x21[i]/d21;
         x32[i] = x32[i]/d32;
         x43[i] = x43[i]/d43;
        }

        ab = dotv(x21,x32);
        bc = dotv(x32,x43);
        ac = dotv(x21,x43);

        vecv(bxc, x32, x43);
        abxc=dotv(x21, bxc);

        t = 180.0 * atan2f(abxc, -ac + ab*bc )/M_PI;
        return t;

}
/*
double torsionv(double *x1, double *x2, double *x3, double *x4)
{
    double d21,d32,d43,x21[3],x32[3],x43[3],bxc[3];
    double  ac,ab,bc,abxc, t;
    int i;

    diffv(x21, x2, x1);
    diffv(x32, x3, x2);
    diffv(x43, x4, x3);
    d21 = sqrt(dotv(x21,x21));
    d32 = sqrt(dotv(x32,x32));
    d43 = sqrt(dotv(x43,x43));

    for(i = 0; i < 3; i++)
    {
     	x21[i] = x21[i]/d21;
        x32[i] = x32[i]/d32;
        x43[i] = x43[i]/d43;
    }

    ab = dotv(x21,x32);
    bc = dotv(x32,x43);
    ac = dotv(x21,x43);

    vecv(bxc, x32, x43);
    abxc = dotv(x21, bxc);

    t = 180.0 * atan2(abxc, -ac + ab*bc )/M_PI;
    return t;
}
*/
/* distanza tra due punti (x1,y1,z1) (x2,y2,z2)*/
double dist(double x1, double y1, double z1, double x2, double y2, double z2)
{
	double d;
	d = dot(x1-x2, y1-y2, z1-z2, x1-x2, y1-y2, z1-z2);
 	d = sqrt(d);
	return d; 
}

/* calcolo angolo (versione coi punti) */
double angle(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3)
{
    double d1,d2,angle;

    d1 = dot(x1-x2, y1-y2, z1-z2,x1-x2, y1-y2, z1-z2);
    d2 = dot(x3-x2, y3-y2, z3-z2,x3-x2, y3-y2, z3-z2);
    angle = dot((x1-x2)/sqrt(d1), (y1-y2)/sqrt(d1), (z1-z2)/sqrt(d1), (x3-x2)/sqrt(d2), (y3-y2)/sqrt(d2), (z3-z2)/sqrt(d2));

    angle = 180.0 * acos(angle) / M_PI;
    return angle;
}

/* calcolo angolo torsionale (versione coi punti)*/
double torsion(double x1, double y1, double z1, 
                            double x2, double y2, double z2,
                            double x3, double y3, double z3, 
                            double x4, double y4, double z4)
{
    double d1,d2,d3,p1_x,p1_y,p1_z,p2_x,p2_y,p2_z,p3_x,p3_y,p3_z;
    double  ac,ab,bc,abxc,bxc_x,bxc_y,bxc_z;
    double torsion;

    d1 = dot(x2-x1,y2-y1,z2-z1, x2-x1,y2-y1,z2-z1);
    d1 =  sqrt(d1);
    d2 = dot(x3-x2,y3-y2,z3-z2, x3-x2,y3-y2,z3-z2);
    d2 =  sqrt(d2);
    d3 = dot(x4-x3,y4-y3,z4-z3, x4-x3,y4-y3,z4-z3);
    d3 =  sqrt(d3);

    p1_x = (x2-x1)/d1; p1_y = (y2-y1)/d1; p1_z = (z2-z1)/d1;
    p2_x = (x3-x2)/d2; p2_y = (y3-y2)/d2; p2_z = (z3-z2)/d2;
    p3_x = (x4-x3)/d3; p3_y = (y4-y3)/d3; p3_z = (z4-z3)/d3;

    ab = dot(p1_x, p1_y, p1_z, p2_x, p2_y, p2_z);
    bc = dot(p2_x, p2_y, p2_z, p3_x, p3_y, p3_z);
    ac = dot(p1_x, p1_y, p1_z, p3_x, p3_y, p3_z);

    vec(&bxc_x, &bxc_y, &bxc_z, p2_x, p2_y, p2_z, p3_x, p3_y, p3_z);
    abxc = dot(p1_x, p1_y, p1_z, bxc_x, bxc_y, bxc_z);

    torsion = 180.0 * atan2f(abxc, -ac + ab*bc )/M_PI;
    return torsion;
}

/* prodotto scalare tra due vettori (versione coi punti)*/
double dot(double x1, double y1, double z1, double x2, double y2, double z2)
{
    double d;
    d = x1*x2 +  y1*y2 + z1*z2;
    return d;
}

/* prodotto vettoriale (versione coi punti)*/
void vec(double *x, double *y, double *z,
         double x1, double y1, double z1, double x2, double y2, double z2)
{
    *x =  (y1 * z2 - z1 * y2);
    *y = -(x1 * z2 - z1 * x2);
    *z =  (x1 * y2 - y1 * x2);
}

/* rotazione di un vettore v attorno ad un asse a per un punto o di un angolo theta */
void rotv(double *o, double *a, double theta, double *v)
{
        double ov[3],x[3],y[3],z[3];
        double  mod,c,s;
        int i;

        mod = sqrt(dotv(a,a));
        for(i=0; i< 3; i++)
         a[i] = a[i]/mod;


        for(i=0; i< 3; i++)
        diffv(ov, v, o); 
        for(i=0; i< 3; i++)
         {
         z[i] = a[i] * dotv(a,ov);
         x[i] = ov[i] - z[i];
         }
        vecv(y,a,x);
        c = cos(M_PI * theta/180.0);
        s = sin(M_PI * theta/180.0);
/*        printf("%11.3f %11.3f %11.3f\n",theta, c, s); */
        for(i=0; i< 3; i++)
         v[i] = x[i] * c + y[i] * s + z[i] + o[i];
}

/*** rotate about the origin and translate a vector ***/
void rot_trans(double *v_in, double **R, double *t, double *v_out)
{
        int i,j;
        double copy_v_in[3];
for(i=0; i<3; i++)
copy_v_in[i] = v_in[i];
for(i=0; i<3; i++)
{
v_out[i] = 0.0;
for(j=0; j<3; j++)
v_out[i] = v_out[i] + R[i][j] * copy_v_in[j];
v_out[i] = v_out[i] + t[i];
}

}


void apply_rot_trans(double **mob, int n, double *t, double **R)
{
/* trova il centro di massa, ruota attorno al centro di massa 
e trasla */
int i, j, k;
double *cm, *v;

cm = (double *) calloc((size_t) 3, sizeof(double));
v = (double *) calloc((size_t) 3, sizeof(double));

for(j=0; j< 3; j++)
for(i=0; i< n; i++)
cm[j] = cm[j] + mob[i][j];

for(j=0; j< 3; j++)
cm[j] = cm[j]/(double) n;

/*
      printf("mob:\n");

      for(i=0; i< n; i++)
      {
      for(j=0; j<3; j++)
      printf("%8.5f ",mob[i][j]);
      printf("\n");
      }

      printf("cm:\n");
      for(j=0; j<3; j++)
      printf("%8.5f ",cm[j]);
      printf("\n");
*/
for(j=0; j< 3; j++)
for(i=0; i< n; i++)
mob[i][j] = -cm[j] + mob[i][j];

/*      printf("mob - cm:\n");
      for(i=0; i< n; i++)
      {
      for(j=0; j<3; j++)
      printf("%8.5f ",mob[i][j]);
      printf("\n");
      }
*/

for(i=0; i< n; i++)
{
for(j=0; j< 3; j++)
{
v[j] = 0;
for(k=0; k< 3; k++)
v[j] = v[j] + R[j][k] * mob[i][k] ;
}
for(j=0; j< 3; j++)
mob[i][j] = v[j] + cm[j] + t[j];
/*      for(j=0; j<3; j++)
      printf("%8.5f ",mob[i][j]);
      printf("\n");
*/

}

free(v);
free(cm);

}

void apply_c_rot_trans(double **mob, int n, double *c, double *t, double **R )
{
int i, j, k;
double *v;

v = (double *) calloc((size_t) 3, sizeof(double));

/*
      printf("mob:\n");

      for(i=0; i< n; i++)
      {
      for(j=0; j<3; j++)
      printf("%8.5f ",mob[i][j]);
      printf("\n");
      }

      printf("cm:\n");
      for(j=0; j<3; j++)
      printf("%8.5f ",cm[j]);
      printf("\n");
*/
for(j=0; j< 3; j++)
for(i=0; i< n; i++)
mob[i][j] = -c[j] + mob[i][j];

/*      printf("mob - cm:\n");
      for(i=0; i< n; i++)
      {
      for(j=0; j<3; j++)
      printf("%8.5f ",mob[i][j]);
      printf("\n");
      }
*/

for(i=0; i< n; i++)
{
for(j=0; j< 3; j++)
{
v[j] = 0;
for(k=0; k< 3; k++)
v[j] = v[j] + R[j][k] * mob[i][k] ;
}
for(j=0; j< 3; j++)
mob[i][j] = v[j] + c[j] + t[j];
/*      for(j=0; j<3; j++)
      printf("%8.5f ",mob[i][j]);
      printf("\n");
*/

}

free(v);

}

void rmsd_no_fit(double **ref, double **mob, int n, double *rmsd)
{
int i,j;
*rmsd = 0.0;
for(i=0; i<n; i++)
for(j=0; j<3; j++)
*rmsd = *rmsd + (ref[i][j] - mob[i][j]) *  (ref[i][j] - mob[i][j]);
*rmsd = sqrt(*rmsd/(double) n);
}

void center_system(struct System *system)
{
double x[3];
x[0] = x[1] = x[2] = 0.0;
int i,j;
for(i=0;i<(*system).n_atoms; i++)
for(j=0;j<3;j++)
x[j] = x[j] + (*system).atoms[i].coor[j];
for(j=0;j<3;j++)
x[j] = x[j] / (double) (*system).n_atoms;
for(i=0;i<(*system).n_atoms; i++)
for(j=0;j<3;j++)
(*system).atoms[i].coor[j] = (*system).atoms[i].coor[j] - x[j];
}

void cmit(struct System system, double *c, double **I)
{
double mass;
int i,j,k;
double v[3];
c[0] = c[1] = c[2] = 0.0;
for(i=0;i<(system).n_atoms; i++)
{
for(j=0;j<3;j++)
c[j] = c[j] + system.atoms[i].coor[j] * system.atoms[i].mass;
mass = mass + system.atoms[i].mass;
}
for(j=0;j<3;j++)
c[j] = c[j] / mass;
for(j=0;j<3;j++)
for(k=0;k<3;k++)
I[j][k] = 0.0;
for(i=0;i<(system).n_atoms; i++)
{
diffv(system.atoms[i].coor, c, v);
for(j=0;j<3;j++)
for(k=0;k<3;k++)
{
I[j][k] = I[j][k] - system.atoms[i].mass * v[j]*v[k];
if(j==k) 
I[j][k] = I[j][k] + system.atoms[i].mass * (v[0]*v[0] + v[1]*v[1] * v[2]*v[2]); 
}
}
}
void translate_system(struct System *system, double *t)
{
int i,j;
for(i=0;i<(*system).n_atoms; i++)
for(j=0;j<3;j++)
(*system).atoms[i].coor[j] = (*system).atoms[i].coor[j] + t[j];
}

/* Scrive sul file fp n_atoms atomi dall'array atoms*/
void write_PDB_atoms(FILE *fp, int n_atoms, struct Atom *atoms, struct Trj trj)
{
	int i,j,k,l;
	char buf[120];

        if(trj.naxf != n_atoms)
          {
          printf("the number of atoms per frame (%i) is not equal to number of atoms (%i)... \n", trj.naxf, n_atoms); 
//          exit(0);
          }
        for(j = 1,l=0; j<= trj.nf; j++)
        {
	fprintf(fp, "MODEL %8i\n", j);
	for(i=0;i< n_atoms; i++)
	{ 
                for(k = 0; k< 3; k++)
                atoms[i].coor[k] = trj.coor[l][k];
		write_atom_pdb(buf, atoms[i]);
		fprintf(fp, "%s", buf);
                l++;
	}
	fprintf(fp, "ENDMDL\n");
        }
}

void write_PQR_atoms(FILE *fp, int n_atoms, struct Atom *atoms, struct Trj trj)
{
	int i,j,k,l;
	char buf[120];

        if(trj.naxf != n_atoms)
          {
          printf("the number of atoms per frame (%i) is not equal to number of atoms (%i)... exiting\n", trj.naxf, n_atoms); 
          exit(0);
          }
        for(j = 1,l=0; j<= trj.nf; j++)
        {
	fprintf(fp, "MODEL %8i\n", j);
	for(i=0;i< n_atoms; i++)
	{ 
                for(k = 0; k< 3; k++)
                atoms[i].coor[k] = trj.coor[l][k];
		write_atom_pqr(buf, atoms[i]);
		fprintf(fp, "%s", buf);
                l++;
	}
	fprintf(fp, "ENDMDL\n");
        }
}

/*Costruisce una stringa buf contenente le informazioni della 
struttura atom secondo la sintassi dei file pdb 
(ATTENZIONE: c'è un parametro di sprintf in più*/
void write_atom_pdb(char *buf, struct Atom atom)
{
	int i;
	char temp[5] = "    ";
	if(strlen(atom.at_name) < 4)
		for(i = 1; i<= strlen(atom.at_name); i++) 
			temp[i] = atom.at_name[i-1]; 
	else strcpy(temp,atom.at_name);
	sprintf(buf,"ATOM  %5i %4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%6.2lf%6.2lf      %4s%2s%2s\n",
             atom.at_n%100000, temp, atom.alt_loc, atom.res_name, atom.chain,
             atom.res_n, atom.res_ins,atom.coor[0], atom.coor[1], atom.coor[2], 
             atom.occ, atom.temp, atom.segid, atom.element, atom.pdb_chrg); 
}

void write_atom_pqr(char *buf, struct Atom atom)
{
	int i;
	char temp[5] = "    ";
	if(strlen(atom.at_name) < 4)
		for(i = 1; i<= strlen(atom.at_name); i++) 
			temp[i] = atom.at_name[i-1]; 
	else strcpy(temp,atom.at_name);
	sprintf(buf,"ATOM  %5i %4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f %7.4lf %7.4lf      %4s%2s%2s\n",
             atom.at_n, temp, atom.alt_loc, atom.res_name, atom.chain,
             atom.res_n, atom.res_ins,atom.coor[0], atom.coor[1], atom.coor[2], 
             atom.charge, atom.radius, atom.segid, atom.element, atom.pdb_chrg); 
}

void copy_atom(struct Atom *atom_to, struct Atom atom_from)
{
int i;
strcpy((*atom_to).at_name,atom_from.at_name);
(*atom_to).alt_loc = atom_from.alt_loc;
strcpy((*atom_to).res_name,atom_from.res_name);
(*atom_to).chain = atom_from.chain;
strcpy((*atom_to).element,atom_from.element);
(*atom_to).model = atom_from.model;
(*atom_to).at_n = atom_from.at_n;
(*atom_to).res_n = atom_from.res_n;
(*atom_to).res_ins = atom_from.res_ins;

for(i=0;i<3;i++)
{
//printf("copio %8.3lf\n", atom_from.coor[i]);
(*atom_to).coor[i] = atom_from.coor[i];
//printf("copiato %8.3lf\n", (*atom_to).coor[i]);
}
(*atom_to).occ = atom_from.occ;
(*atom_to).temp = atom_from.temp;
strcpy((*atom_to).segid,atom_from.segid);
strcpy((*atom_to).pdb_chrg,atom_from.pdb_chrg);
(*atom_to).mass = atom_from.mass;
strcpy((*atom_to).chem_group,atom_from.chem_group);
strcpy((*atom_to).group,atom_from.group);
(*atom_to).radius = atom_from.radius;
(*atom_to).charge = atom_from.charge;
(*atom_to).atom_type = atom_from.atom_type;
(*atom_to).i = atom_from.i;
}

void read_PQR_atoms(FILE *fp1, int *n_atoms, struct Atom *(*atoms), struct Trj *trj)
{
	char buf[120];
	int i=0, k, n_models,is_trj;
        int mod_id=0;
        struct Atom tmp_atom;

// first check if it is a trajectory or a single structure
	while(fgets(buf,120,fp1) != NULL )
        {
    	if(!strncmp("ATOM",buf,4) && mod_id==0) i++;
	      if(!strncmp("ENDMDL",buf,6)) 
              mod_id++; 
        }
*n_atoms = i; 
n_models = mod_id;
if(mod_id <= 1) mod_id = 1;
  (*atoms) = (struct Atom *) calloc((size_t) *n_atoms , sizeof(struct Atom));
  if((*atoms) == NULL) 
    {
     printf("could not allocate memory for %i atoms... exiting...\n", *n_atoms);
     exit(0);
    }
  (*trj).nf = mod_id; 
  (*trj).naxf = *n_atoms; 
  trj->coor = (double **) calloc((size_t) ((*trj).nf * (*trj).naxf), sizeof(double *));
  if(trj->coor == NULL) 
    {
     printf("could not allocate memory for %i atoms... exiting...\n", (*trj).nf * (*trj).naxf);
     exit(0);
    }
   for(i = 0; i < (*trj).nf * (*trj).naxf; i++)
    {
    trj->coor[i] =  (double *) calloc((size_t) 3, sizeof(double));
    if(trj->coor[i] == NULL) 
    {
     printf("could not allocate memory for %i-th atom coordinates... exiting...\n", i);
     exit(0);
    }
    }
        rewind(fp1);
        i = 0;
	while(fgets(buf,120,fp1) != NULL)
    	if(!strncmp("ATOM",buf,4)) 
            {
                        if(i < *n_atoms)
                        {
			read_atom_pqr(buf, &((*atoms)[i]));
                        for(k=0;k<3;k++)
                        (*trj).coor[i][k] = (*atoms)[i].coor[k];
                        i++;
                        }
                        else if(mod_id > 1)
                        {
			read_atom_pqr(buf, &tmp_atom);
                        for(k=0;k<3;k++)
                        (*trj).coor[i][k] = tmp_atom.coor[k];
                        i++;
                        }
            if(!(i%10000))  printf("%i atoms read\n",i);
	    }
	    else
	      if(!strncmp("ENDMDL",buf,6)) 
              mod_id++; 
}

/* Legge un le informazioni di un atomo da una stringa buf e le salva 
nella strutture atom.*/
void read_atom_pqr(char *buf, struct Atom *atom)
{

	char at_rec[5];
    char tok[10];              
 
    strncpy(tok,buf,4);
    tok[4] = '\0';
    sscanf(tok,"%s", at_rec);
    if(strncmp("ATOM",at_rec,4))
    {
    	printf("The ATOM line does not start with string ATOM... exiting...\n");
        exit(1);
    }

    strncpy(tok,buf + 6,5);
    tok[5] = '\0';
    sscanf(tok,"%i",&(atom->at_n));

    strncpy(tok,buf + 12,4);
    tok[4] = '\0';
    sscanf(tok,"%s", atom->at_name);
 
    strncpy(tok,buf + 16,1);
    tok[1] = '\0';
    if(sscanf(tok,"%c", &(atom->alt_loc)) == -1) atom->alt_loc=' ';
/*    else if ((atom->alt_loc=='A') || (atom->alt_loc=='1')) 
    atom->alt_loc=' '; */ 

	strncpy(tok,buf + 17,3);
    tok[3] = '\0'; 
    sscanf(tok,"%s", atom->res_name);

	strncpy(tok,buf + 21,1);
    tok[1] = '\0';
    if(sscanf(tok,"%c", &(atom->chain)) == EOF) atom->chain = ' ';

    strncpy(tok,buf + 22,4);
    tok[4] = '\0';
    sscanf(tok,"%i", &(atom->res_n));

	strncpy(tok,buf + 26,1);
    tok[1] = '\0';
    if (sscanf(tok,"%c", &(atom->res_ins)) == EOF) atom->res_ins=' ';

    strncpy(tok,buf + 30,8);
    tok[8] = '\0';
    sscanf(tok,"%lf", &(atom->coor[0]));

	strncpy(tok,buf + 38,8);
    tok[8] = '\0';
    sscanf(tok,"%lf", &(atom->coor[1]));

    strncpy(tok,buf + 46,8);
    tok[8] = '\0';
    sscanf(tok,"%lf", &(atom->coor[2]));

    sscanf(buf + 55 ,"%lf %lf", &(atom->charge), &(atom->radius));

/* the following is a patch just to read pqr from pdb2pqr.py */


/*
    strncpy(tok,buf + 54,6);
    tok[6] = '\0';
    sscanf(tok,"%lf", &(atom->occ));

	strncpy(tok,buf + 60,6);
    tok[6] = '\0';
    sscanf(tok,"%lf", &(atom->temp));

    if(strlen(buf) > 72)
    { 
      	strncpy(tok,buf + 72,4);
       	tok[4] = '\0';
		if (sscanf(tok,"%s", (atom->segid)) == EOF) 
			strcpy(atom->segid,"    ");
    }
    else strcpy(atom->segid,"    ");


    if(strlen(buf) > 76)
	{
       	strncpy(tok,buf + 76,2);
       	tok[2] = '\0';
       	if (sscanf(tok,"%s", (atom->element)) == EOF) 
			strcpy(atom->element,"UN");
	}

    if(strlen(buf) > 78)
	{
       	strncpy(tok,buf + 78,2);
       	tok[2] = '\0';
       	if (sscanf(tok,"%s", (atom->pdb_chrg)) == EOF) 
			strcpy(atom->pdb_chrg,"  ");
	}
*/
}


struct Neighbour *neighbour_alloc(int n, int max_nei)
{
    int i;
    struct Neighbour *p;
    p = (struct Neighbour *) calloc((size_t) n , sizeof(struct Neighbour));
    if(p==NULL) printf("Could not allocate memory for Neighbour... Exiting...\n");
//    printf("%i %i\n", n, max_nei);
    for(i = 0; i< n; i++)
    p[i].max_n = max_nei;

    for(i = 0; i< n; i++)
    {
    p[i].list = (int *) calloc((size_t) max_nei, sizeof(int));
    if(p[i].list == NULL) 
    { printf("Could not allocate memory Neighbour.list[%i]... Exiting...\n", i);    exit(0);}
    }
    for(i = 0; i< n; i++)
    {
    p[i].d = (double *) calloc((size_t) max_nei, sizeof(double));
    if(p[i].d == NULL) 
    { printf("Could not allocate memory Neighbour.d[%i]... Exiting...\n", i);    exit(0);}
    }
    return p;
}

void neighbour_free(struct Neighbour *neighbours, int n)
{
    int i;
    for(i = 0; i< n; i++)
    free((neighbours)[i].list);
    free(neighbours);
}

void neighbour_realloc(struct Neighbour *neighbours, int n)
{
    int i;
    for(i = 0; i< n; i++)
    {
    neighbours[i].max_n = neighbours[i].n_neighbours;
    neighbours[i].list = (int *) realloc(neighbours[i].list, (size_t) neighbours[i].n_neighbours * sizeof(int));
    neighbours[i].d = (double *) realloc(neighbours[i].d, (size_t) neighbours[i].n_neighbours * sizeof(double));
    }
}


void atom_grid_free(struct Atom_grid *atom_grid)
{
int i;
for(i=0; i< (*atom_grid).grid_size; i++)
    free((*atom_grid).atom_node[i]);
free((*atom_grid).atom_node); 
free((*atom_grid).n_atom_node); 
(*atom_grid).max_atom_node = 0;
(*atom_grid).grid_size = 0;	
}

void atom_grid_alloc(struct Atom_grid *atom_grid)
{
int i;
//printf("Grid size: %i\n", (*atom_grid).grid_size);
atom_grid->atom_node = (int **) calloc((size_t) (*atom_grid).grid_size , sizeof(int*));
if(atom_grid->atom_node == NULL)
    { printf("Could not allocate memory for atom_grid->atom_node... Exiting...\n");    exit(0);}
//printf("d(n_atom_node) = %i\n", (*atom_grid).grid_size); 
atom_grid->n_atom_node = (int *) calloc((size_t) (*atom_grid).grid_size , sizeof (int ));
if(atom_grid->n_atom_node == NULL)
    { printf("Could not allocate memory for atom_grid->n_atom_node... Exiting...\n");    exit(0);}
for(i = 0; i< (*atom_grid).grid_size; i++)
 {
(*atom_grid).atom_node[i] = (int *) calloc((size_t) (*atom_grid).max_atom_node , sizeof(int));
if(atom_grid->atom_node[i] == NULL)
    { printf("Could not allocate memory for atom_grid->atom_node[%i]... Exiting...\n", i);    exit(0);}
 }
atom_grid->PBC_ON = 0;
}

void atom_grid_realloc(struct Atom_grid *atom_grid)
{
int i;
for(i=0; i< (*atom_grid).grid_size;i++)
atom_grid->atom_node[i] = (int *) realloc(atom_grid->atom_node[i], (size_t) (*atom_grid).n_atom_node[i] * sizeof(int));
}

FILE *file_open(char *fname,char *acc) {
    FILE *fp;
    fp =fopen(fname,acc);
    if (fp==NULL) 
	{
        fprintf(stderr,"unable to open file %s\n",fname);
        exit(1);
    }
    return(fp);
}

/************** external variables *********************/

/***** local includes *********************/

//#include <pdblib.h>
/**************************************************
*  Each program has its options in flag_par 
*  (flags or values for parameters)
*
*  init_flag_par initializes flag_par
*  check_cmd_line rears options from command line
*  print_info_flag_par prints options       
*
***************************************************/
struct Flag_par {  
	       char file_in_pdb[256];
	       char file_out[256];
	       char file_msms[256];
               double msms_area;
               double sasa_area;
               double probe_radius;
               double test_charge_radius;
               double min_radius;
               double salt_radius;
               double pdie;
               double sdie;
               double ions;
               double temp;
               double cutoff;
               double cutoff2;
               double maxf;
               double gamma;
               double kp;
               int gbr;
               int surface;
               int srfpot;
               int srfpotav;
               int solv_e;
               int pqg;
               int grid, nx, ny, nz;
	       double mesh;
               int msms;
               int pka;
               int autoter;
               char pkadef[256];
               int defpka;
               double pkadie;
               int pkass;
               int verbose;
	     } Flag_par;

struct Area_atom {
        double ses, sas;
        double area;
        int n_model;
        char  segment[5]; 
        char chain; 
	int res_n;
	char res_ins; 
        char  at_name[5]; 
        char alt_loc;
        char res_name[4];
        int at_n;
        } Area_atom;

void grid_srf_dens(
                struct System system,
                struct Atom_grid atom_grid,
                struct Flag_par flag_par,
                struct Srf *srf);

struct Msms_face
{
int n1,n2,n3,type;
};

struct Energy
{
double area, solv, coul, pol, born, *solv_a;
};

struct Const
{
double q, N_av, e0, kb, RT, kd, k_el;
};

struct pka_def
{
int state;
double pka,charge,charge_add;
char atom[5];
char res[5];
};

void get_srf(struct System sys, struct Atom_grid *patom_grid, struct Srf *psrf_ses, struct Srf *psrf_sas, struct Flag_par flag_par);
//void get_srf(struct System sys, struct Atom_grid *patom_grid, struct Srf *psrf_ses, struct At_srf **pses_srf_at, struct Srf *psrf_sas, struct At_srf **psas_srf_at, struct Flag_par flag_par);
	
void pqr2gbr6(struct System sys, struct Srf srf_ses, struct Srf Srf_sas, struct Flag_par flag_par, struct Neighbour *neighbours, double *gbr6_ses, double *gbr6_sas);

void print_radii(struct System sys, struct Trj trj, double *gbr6_ses, double *gbr6_sas, struct Flag_par flag_par);

void print_srf_pot(struct System sys, struct Trj trj, double *gbr6, struct Srf srf_ses, struct Srf srf_sas, struct Neighbour *neighbours, struct Area_atom *area_atom, struct Flag_par flag_par);

void srfpot(struct System *sys, struct Srf srf, double *gbr6, struct Neighbour *neighbours, struct Flag_par flag_par);

void calc_nrg(struct System sys, struct Energy *pnrg, double *gbr6, struct Srf srf, struct Neighbour *neighbours, struct Const constants, struct Flag_par flag_par);

void calc_print_pka(struct System sys, struct Trj trj, double *gbr6, struct Neighbour *neighbours, struct Const constants, struct Flag_par *flag_par);

void print_nrg(struct Energy nrg, struct System sys, struct Flag_par flag_par);

void out_grid(struct System sys, struct Atom_grid atom_grid, double *gbr6, struct Neighbour *neighbours, struct Const constants, struct Flag_par flag_par);

void scalev(double t, double *v);

double fgb(double rgb1, double rgb2, double rw, double rs, double d, 
           double pdie, double sdie, double kd, double kp);

void mc(double *state, double *site_charge, double *pka, double *de, double **g, double **cova, int n_sites, struct Neighbour *nei_sites, double T, double pH, int nsteps, int presteps, double *energy, double *r_acc);

int cmp_area_atom_by_string(const void *p1, const void *p2);
int cmp_area_atom_by_n(const void *p1, const void *p2);
int cmp_int(const void *p1, const void *p2);

void init_flag_par(struct Flag_par *flag_par);
void check_cmd_line(int argc, char *argv[], struct Flag_par *flag_par);
void print_info_flag_par(struct Flag_par flag_par);

//struct Neighbour ** nei_realloc(struct Neighbour **neighbour_p, int n);
/**** MAIN *****/

int main(int argc, char *argv[]) {
      FILE *fp;
      struct Flag_par flag_par;
      struct Const constants;
      struct System sys;
      struct Trj trj;
      struct Atom_grid atom_grid;
      struct Neighbour *neighbours; 
      struct Energy nrg;
      int i,k; 
      double *gbr6_ses, *gbr6_sas;
      struct Srf *psrf_ses, *psrf_sas;
      struct Srf srf_ses;
      struct Srf srf_sas;
      double charge; 
      struct Area_atom *area_atom;

        /* Initialize options and parameters for this program */
        init_flag_par(&flag_par);      
        /* Read options from command line */
        check_cmd_line(argc, argv, &flag_par);      
	if(flag_par.pka)
	{
		printf("pka calculation selected. Resetting SASA area per point to msms area per point\n");
		printf("this is due to the way pka is calculated....\n");
                flag_par.sasa_area = flag_par.msms_area;
	}
        /* Print information on options and parameters */
        print_info_flag_par(flag_par);      

      constants.q = 1.602176487E-19;
      constants.N_av = 6.02214179E23;
      constants.e0 = 8.8541878176E-12;
      constants.kb = 1.3806504E-23;
      constants.RT = 1.3806504 * 6.02214179 * flag_par.temp;
      constants.kd = 1E-10 * sqrt((1000 * 2 * flag_par.ions * constants.N_av * constants.q * constants.q )/ 
                         (flag_par.sdie * constants.e0 * constants.kb * flag_par.temp));
      constants.k_el = constants.N_av * 1E10 * constants.q * constants.q /(4 * M_PI * constants.e0);
      printf("k_el = %e (kJ/mol) %e (kcal/mol)\n", 
      constants.k_el/1000, constants.k_el/4184); 
      printf("Debye length = %e (A)\n", 1.0/constants.kd);

        /* read the pqr file */
        fp = file_open(flag_par.file_in_pdb,"r");
        read_PQR_atoms(fp, &(sys.n_atoms), &(sys.atoms), &trj);
        fclose(fp);
        printf("\n%i total atoms read\n", sys.n_atoms);

        for(i = 0; i < sys.n_atoms; i++)
        for(k = 0; k < 3; k++)
        sys.atoms[i].coor[k]  = trj.coor[i][k];
        /* reset radius to minimum radius */
        for(i=0; i< sys.n_atoms; i++)
        if(sys.atoms[i].radius < flag_par.min_radius) sys.atoms[i].radius = flag_par.min_radius;

        if(flag_par.verbose)
        for(i=0; i < sys.n_atoms; i++) 
        printf("atom %i: charge %lf  radius %lf\n", i, sys.atoms[i].charge, sys.atoms[i].radius);
// Alloca subito le strutture necessarie
get_grid_parameters(sys, flag_par.probe_radius, flag_par.cutoff, &atom_grid, 'v');
//print_grid_info(atom_grid);
//printf("%i %lf %lf %i %i %i %i %lf\n", sys.n_atoms, flag_par.probe_radius, flag_par.cutoff, atom_grid.grid_X, atom_grid.grid_Y, atom_grid.grid_Z, atom_grid.grid_size, atom_grid.mesh);
atom_grid.PBC_ON = 0;
atom_grid_alloc(&atom_grid);
printf("Initial grid allocation for %i atoms: %i nodes\n", sys.n_atoms, atom_grid.max_atom_node);
neighbours = neighbour_alloc(sys.n_atoms, atom_grid.max_atom_node);
        if(flag_par.srfpot || flag_par.surface)
        area_atom=(struct Area_atom *) calloc((size_t) sys.n_atoms, sizeof(struct Area_atom));

        if(flag_par.solv_e)
        nrg.solv_a=(double *) calloc((size_t) sys.n_atoms, sizeof(double));
        
        gbr6_ses = (double *) calloc((size_t) sys.n_atoms, sizeof(double));
        gbr6_sas = (double *) calloc((size_t) sys.n_atoms, sizeof(double));


	/* print info on total charge */
          charge=0.0;
          for(i=0; i< sys.n_atoms; i++)
               charge = charge + sys.atoms[i].charge;
printf("Total charge on molecule: %11.3f\n", charge);
printf("--------------------------------\n");

//for(i = 0; i <  sys.n_atoms ; i++)
//printf("%5i %i\n", i, srf.at_srf[i].n_isrf);
atoms_on_grid(&atom_grid, sys);
atom_grid_realloc(&atom_grid); 
grid_neighbour(atom_grid, sys, flag_par.cutoff, neighbours,'v');
neighbour_realloc(neighbours, sys.n_atoms);
psrf_ses= &srf_ses;
psrf_sas= &srf_sas;
get_srf(sys, &atom_grid, psrf_ses, psrf_sas, flag_par);
// if msms is not used gbr from SES are approximated
pqr2gbr6(sys, srf_ses, srf_sas, flag_par, neighbours, gbr6_ses, gbr6_sas);
print_radii(sys, trj, gbr6_ses, gbr6_sas, flag_par);
srfpot(&sys, srf_sas, gbr6_ses, neighbours, flag_par); // compute the potential at SAS
print_srf_pot(sys, trj, gbr6_ses, srf_ses, srf_sas, neighbours, area_atom, flag_par);
calc_nrg(sys,&nrg, gbr6_ses, srf_sas, neighbours, constants, flag_par);
print_nrg(nrg,sys,flag_par);
out_grid(sys,atom_grid,gbr6_ses, neighbours, constants,flag_par);
calc_print_pka(sys, trj, gbr6_sas, neighbours, constants, &flag_par);
}

void check_cmd_line(int argc, char *argv[], struct Flag_par *flag_par)
{
	int i;
        char tmp[100];
        char extension[100];

	if(argc < 3) 
	{
	printf("Usage:\n"); 
	printf("bluues2 filename.pqr basename [Options]\n"); 
	printf("Options:\n"); 
	printf("-v (verbose mode)\n"); 
	printf("-m filename.vert (msms .vert file with 3 lines header)\n"); 
	printf("-pa x (msms area per point on SES, 0.1 default)\n"); 
	printf("-p x (use a probe radius for surfacing of x, 1.5 A default)\n"); 
	printf("-tcr x (use a test charge radius of x, 0.7 A default)\n"); 
	printf("-mr x (minimum radius for surfacing of x, 1.0 A default)\n"); 
	printf("-s x (use a salt radius of x, 2.0 default)\n"); 
	printf("-pd x (inner dielectric constant, 4.0 default)\n"); 
	printf("-sd x (outer dielectric constant, 78.54 default)\n"); 
	printf("-kp x (factor for Still formula, 4 default)\n"); 
	printf("-i x (ionic strength (M), 0.15 M default)\n"); 
	printf("-t x (temperature, 298.15 default)\n"); 
	printf("-c x (cutoff for GBR6 calculation), 20.0 A default\n"); 
	printf("-c2 x (cutoff short for GBR6 calculation, 8.0 A default)\n"); 
	printf("-g x (surface tension coefficient, 0.12 J/(A^2 mol) default)\n"); 
	printf("-pka (compute pkas)\n"); 
	printf("-pkadef filename (file containing pKa definitions)\n"); 
	printf("-ter (do not consider terminal N/C as titratable sites)\n"); 
	printf("-pkd x (inner dielectric constant for pka calculations, 20.0 default)\n"); 
	printf("-pss x (number of MC steps per sites in pka calculation, 1000 default)\n"); 
	printf("-dx (output potential grid in dx format)\n"); 
	printf("-nx x (number of grid points in x, 97 default)\n"); 
	printf("-ny x (number of grid points in y, 97 default)\n"); 
	printf("-nz x (number of grid points in z, 97 default)\n"); 
	printf("-mesh x (grid mesh size, 1.0 default)\n"); 
	printf("-srf (output surface points, potential, atom n. and normal vectors)\n"); 
	printf("-srfpot (output average atomic surface potential)\n"); 
	printf("-srfpotav (output atomic surface potential over neighbours within c2 A)\n"); 
	printf("\n"); 
	exit(1);
	}
        
        strcpy((*flag_par).file_in_pdb, argv[1]);
        strcpy((*flag_par).file_out, argv[2]);

	for (i = 3; i < argc; i++) 
        {
		if (!strncmp(argv[i],"-srfpot",8)) (*flag_par).srfpot = 1; 
		else if (!strcmp(argv[i],"-srfpotav")) 
                    {
                      (*flag_par).srfpot = 1; 
                      (*flag_par).srfpotav = 1;
		    } 
		else if (!strncmp(argv[i],"-pka",5)) (*flag_par).pka = 1; 
		else if (!strncmp(argv[i],"-ter",5)) (*flag_par).autoter = 0; 
		else if (!strncmp(argv[i],"-pkd",5)) sscanf(argv[++i], "%lf", &((*flag_par).pkadie)); 
		else if (!strcmp(argv[i],"-pkadef")) {sscanf(argv[++i], "%s", (*flag_par).pkadef); (*flag_par).defpka = 1;}
		else if (!strncmp(argv[i],"-srf",5)) (*flag_par).surface = 1; 
		else if (!strncmp(argv[i],"-pss",5)) sscanf(argv[++i], "%i", &((*flag_par).pkass)); 
		else if (!strncmp(argv[i],"-tcr",5)) sscanf(argv[++i], "%lf", &((*flag_par).test_charge_radius)); 
		else if (!strncmp(argv[i],"-pd",4)) sscanf(argv[++i], "%lf", &((*flag_par).pdie)); 
		else if (!strncmp(argv[i],"-sd",4)) sscanf(argv[++i], "%lf", &((*flag_par).sdie)); 
		else if (!strncmp(argv[i],"-mr",4)) sscanf(argv[++i], "%lf", &((*flag_par).min_radius)); 
                else if (!strcmp(argv[i],"-pases")) sscanf(argv[++i], "%lf", &((*flag_par).msms_area));
                else if (!strcmp(argv[i],"-pasas")) sscanf(argv[++i], "%lf", &((*flag_par).sasa_area));
                else if (!strncmp(argv[i],"-dx",4)) (*flag_par).grid = 1;
		else if (!strncmp(argv[i],"-mesh",6)) sscanf(argv[++i], "%lf", &((*flag_par).mesh)); 
		else if (!strncmp(argv[i],"-nx",4)) sscanf(argv[++i], "%i", &((*flag_par).nx)); 
		else if (!strncmp(argv[i],"-ny",4)) sscanf(argv[++i], "%i", &((*flag_par).ny)); 
		else if (!strncmp(argv[i],"-nz",4)) sscanf(argv[++i], "%i", &((*flag_par).nz)); 
		else if (!strncmp(argv[i],"-i",3)) sscanf(argv[++i], "%lf", &((*flag_par).ions)); 
		else if (!strncmp(argv[i],"-kp",4)) sscanf(argv[++i], "%lf", &((*flag_par).kp)); 
		else if (!strncmp(argv[i],"-c2",4)) sscanf(argv[++i], "%lf", &((*flag_par).cutoff2)); 
 		else if (!strncmp(argv[i],"-m",3)) {sscanf(argv[++i], "%s", ((*flag_par).file_msms)); (*flag_par).msms = 1;}
		else if (!strncmp(argv[i],"-p",3)) sscanf(argv[++i], "%lf", &((*flag_par).probe_radius)); 
		else if (!strncmp(argv[i],"-s",3)) sscanf(argv[++i], "%lf", &((*flag_par).salt_radius)); 
		else if (!strncmp(argv[i],"-t",3)) sscanf(argv[++i], "%lf", &((*flag_par).temp)); 
		else if (!strncmp(argv[i],"-g",3)) sscanf(argv[++i], "%lf", &((*flag_par).gamma)); 
		else if (!strncmp(argv[i],"-i",3)) sscanf(argv[++i], "%lf", &((*flag_par).ions)); 
		else if (!strncmp(argv[i],"-c",3)) sscanf(argv[++i], "%lf", &((*flag_par).cutoff)); 
		else if (!strncmp(argv[i],"-g",3)) (*flag_par).grid = 1; 
		else if (!strncmp(argv[i],"-v",3)) (*flag_par).verbose = 1;
		else 
                {
                 printf("I don't know option %s\n", argv[i]);
                 exit(2);
                }
        }

}

void init_flag_par(struct Flag_par *flag_par)
{
(*flag_par).verbose=0;
(*flag_par).probe_radius=1.5;
(*flag_par).test_charge_radius=0.7;
(*flag_par).min_radius=1.0;
(*flag_par).salt_radius=2.0;
(*flag_par).ions=0.150;
(*flag_par).pdie=4.0;
(*flag_par).sdie=78.54;
(*flag_par).kp=4;
(*flag_par).temp=298.15;
(*flag_par).cutoff=20.0;
(*flag_par).cutoff2=8.0;
(*flag_par).gbr = 1;
(*flag_par).srfpot = 0;
(*flag_par).srfpotav = 0;
(*flag_par).solv_e = 1;
(*flag_par).maxf = 999.94999;
(*flag_par).gamma = 0.12;
(*flag_par).grid = 0;
(*flag_par).pka = 0;
(*flag_par).autoter = 1;
(*flag_par).pkass = 1000;
(*flag_par).pkadie=20.0;
(*flag_par).defpka=0;
(*flag_par).msms=0;
(*flag_par).msms_area=0.1;
(*flag_par).sasa_area=0.2;
(*flag_par).surface=0.0;
(*flag_par).pqg = 1;
(*flag_par).nx = 97;
(*flag_par).ny = 97;
(*flag_par).nz = 97;
(*flag_par).mesh = 1.0;
}

void print_info_flag_par(struct Flag_par flag_par)
{
        printf("--------------------------------\nRUN PARAMETERS:\n--------------------------------\n");
        printf("pdb file: %s\n", flag_par.file_in_pdb);
        printf("output file: %s\n", flag_par.file_out);
        printf("probe_radius: %lf\n", flag_par.probe_radius);
        printf("minimum atomic radius: %lf\n", flag_par.min_radius);
        printf("ionic strength (M): %lf\n", flag_par.ions);
        printf("salt_radius (A): %lf\n", flag_par.salt_radius);
        printf("inner dielectric constant: %lf\n", flag_par.pdie);
        printf("outer dielectric constant: %lf\n", flag_par.sdie);
        printf("temperature (K): %lf\n", flag_par.temp);
        printf("cutoff for GBR6 calculation (A): %lf\n", flag_par.cutoff);
        printf("cutoff for electrostatic patches (A): %lf\n", flag_par.cutoff2);
        printf("gamma (kJ/A^2): %lf\n", flag_par.gamma);
        if(flag_par.pka)
        {
        printf("inner dielectric constant in pka calculations: %lf\n", flag_par.pkadie);
        printf("pka steps per site: %i\n", flag_par.pkass);
        if(flag_par.defpka)
		printf("I will read pKa definitions in %s\n", flag_par.pkadef);
       if(flag_par.autoter)
		printf("I will consider N-ter N and C-ter C as ionizable termini\n");
        }
        if(flag_par.srfpot || flag_par.grid)
        {
        printf("test charge radius: %lf\n", flag_par.test_charge_radius);
        }
        if(flag_par.msms )
        {
        printf("input file from MSMS: %s\n", flag_par.file_msms);
        }
        printf("area (A^2) per point on SES: %lf\n", flag_par.msms_area);
        printf("--------------------------------\n");
}

double fgb(double rgb1, double rgb2, double rw, double rs, double d, 
           double pdie, double sdie, double kd, double kp)
{
double f, rs_eff, rgb, deff;
rgb = sqrt(rgb1 * rgb2);
if(rs < rw) rs = rw; // aggiunto per evitare nan
rs = sqrt(
(rgb1 + rs - rw) *
(rgb2 + rs - rw)
);
deff = sqrt(d*d + rgb * rgb *exp(-d*d/(kp * rgb * rgb)));
//printf("DEBUG: rgb = %e   rgb 1 = %e    rgb 2 = %e   d_eff = %e k_d = %e rs = %e \n", rgb, rgb1, rgb2, deff, kd, rs);
//printf("DEBUG: exp = %e   den = %e\n", (exp(-kd * (deff - rs))), (deff * (1 + kd * rs)));
if(deff >= rs)
f = -1/(pdie * deff) + (1/sdie) * (exp(-kd * (deff - rs)) /(deff * (1 + kd * rs)));
else 
f = -1/(pdie * deff) + (1/sdie) * (1.0/(rs * (1 + kd * rs)) + 1.0/deff - 1.0/rs);
/*printf("rgb1 = %lf rgb2 = %lf dist = %lf deff = %lf deff_attesa = %lf\n", 
rgb1, rgb2, d, deff, sqrt(rgb1*rgb2 + d *d)); */
return f;
}

void mc(double *state, double *site_charge, double *pka, double *dpka, double **g, double **cova, int n_sites, struct Neighbour *nei_sites, double T, double pH, int nsteps, int presteps, double *energy, double *r_acc)
{
int i,j,k,l,ll,n,n_accept,n_max;
double f, f_new, tmpf, lRT, RT;
double *a_new, *a_ave, f_ave, f_0;
n_max = n_sites*5;
lRT = log(10) * 6.02214179 * 1.3806504 * T;
RT = 6.02214179 * 1.3806504 * T;
a_new = (double *) calloc((size_t) n_sites, sizeof(double));
a_ave = (double *) calloc((size_t) n_sites, sizeof(double));
f_ave = 0.0;

for(i=0; i<n_sites; i++)
if(state[i] > 0.5) state[i] = 1.0;
else state[i] = 0.0;

for(k=0; k<n_sites; k++) a_new[k] = state[k];
f = 0;
for(k=0; k<n_sites; k++)
{
f = f + state[k] * site_charge[k] * (+pH - (pka[k] + dpka[k])) * lRT;
for(ll=0; ll<nei_sites[k].n_neighbours; ll++) 
{
l = nei_sites[k].list[ll];
if(l>k)
f = f + state[k] * state[l] * g[k][l]; 
}
}
n_accept = 0;
for(k=0; k<presteps; k++)
{
l = 0;
do
{
i = rand()%n_sites;
//printf("pH = %f pka = %f %i mosse\n", pH, pka[i], l);
l++;
}
while(((fabs(pka[i] - pH)  > 6.0) && l < n_max) && l< 1);
//Mossa
tmpf = state[i];
if(state[i] > 0.5) state[i] = 0.0;
else state[i] = 1.0;
f_new = f + (state[i] - tmpf)* ( site_charge[i] * (+pH - (pka[i] + dpka[i])) * lRT);
for(ll=0; ll<nei_sites[i].n_neighbours; ll++) 
{
l = nei_sites[i].list[ll];
if(l!=i)
f_new = f_new + (state[i] - tmpf) * state[l] * g[i][l]; 
}
if((f_new < f) || exp( -(f_new-f) /RT ) > drand48())
f = f_new;
else
state[i] = tmpf;
}
for(k=0; k<nsteps; k++)
{
l = 0;
do
{
i = rand()%n_sites;
l++;
}
while(((fabs(pka[i] - pH)  > 6.0) && l < n_max) && l < 1);
tmpf = state[i];
if(state[i] > 0.5) state[i] = 0.0;
else state[i] = 1.0;
f_new = f + (state[i] - tmpf) * (site_charge[i] * (-(pka[i] + dpka[i]) + pH) * lRT);

for(ll=0; ll<nei_sites[i].n_neighbours; ll++) 
{
l = nei_sites[i].list[ll];
if(l!=i)
f_new = f_new + (state[i] - tmpf) * state[l] * g[i][l]; 
}

if((f_new < f) || exp( -(f_new-f) /RT ) > drand48())
{
f = f_new;
n_accept++;
}
else
state[i] = tmpf;
//Medie
for(j=0; j< n_sites; j++)
{
a_ave[j] = a_ave[j] + state[j];
//for(l=0; l< n_sites; l++)
//cova[j][l] = cova[j][l] + state[j]*state[l];
}
f_ave = f_ave + f;
}
f_0 = 0.0;
f_ave = f_ave/(double) nsteps;
for(j=0; j< n_sites; j++)
{
if(site_charge[j] > 0.5)
{
tmpf = (1/(pow(10, pH - pka[j]) + 1.0));
f_0 = f_0 + tmpf * ((pH -(pka[j])) * lRT);
}
else
{
tmpf = (1/(pow(10, -pH + pka[j]) + 1.0));
f_0 = f_0 + tmpf * (site_charge[j] * (pH -(pka[j])) * lRT);
}
//printf("tmpf %e  f0 = %e f_ave = %e\n", tmpf, f_0/4186.7999409, f_ave/4186.7999409);
}
*energy = (f_ave - f_0)/1000.0;
//printf("f0 = %e f_ave = %e\n", f_0/4186.7999409, f_ave/4186.7999409);

//Operazioni finali
for(j=0; j< n_sites; j++)
state[j] = a_ave[j]/(double) nsteps;

//for(j=0; j< n_sites; j++)
//for(ll=0; ll<nei_sites[j].n_neighbours; ll++) 
//{
//l = nei_sites[j].list[ll];
//cova[j][l] = cova[j][l]/(double) nsteps  - state[j] * state[l];
//}

*r_acc =  (double) n_accept/ (double) nsteps;
free(a_new);
free(a_ave);
}

int cmp_int(const void *p1, const void *p2)
{
int i1, i2;
i1 = *((int *)p1);
i2 = *((int *)p2);
if(i2>i1) return -1;
else if(i1==i2) return 0;
else return 1;
}

int cmp_area_atom_by_string(const void *p1, const void *p2)
{
	struct Area_atom A, B;
	int check = 0 ; 

	A = *((struct Area_atom *)p1); 
	B = *((struct Area_atom *)p2); 

	if( A.n_model < B.n_model) check = -1; 
        else if ( strcmp(A.segment,B.segment) < 0) check = -1;
        else if ( A.chain < B.chain) check = -1;
        else if ( A.res_n < B.res_n) check = -1;
        else if ( A.res_ins < B.res_ins) check = -1;
        else if ( A.alt_loc < B.alt_loc) check = -1;
        else check = 1;
        return check;
}

int cmp_area_atom_by_n(const void *p1, const void *p2)
{
	struct Area_atom A, B;
	int check = 0 ; 

	A = *((struct Area_atom *)p1); 
	B = *((struct Area_atom *)p2); 

	if( A.at_n < B.at_n) check = -1; 
        else check = 1;
        return check;
}

void pqr2gbr6(struct System sys, struct Srf srf_ses, struct Srf srf_sas, struct Flag_par flag_par, struct Neighbour *neighbours, double *gbr6_ses, double *gbr6_sas)
{
int i, j, k, l;
double *phi,*ir1,*ir2,*ir3, dx = 1e-3;
double *v, tmpf, t0,t1,t2,t3, min[3], max[3], rave, *gbr6;
struct Srf *srf;

        v = (double *) calloc((size_t) 3,sizeof(double));
        phi=(double *) calloc((size_t) sys.n_atoms, sizeof(double));
        ir1=(double *) calloc((size_t) sys.n_atoms, sizeof(double));
        ir2=(double *) calloc((size_t) sys.n_atoms, sizeof(double));
        ir3=(double *) calloc((size_t) sys.n_atoms, sizeof(double));
for(k=0; k<3; k++)
{
min[k] =  99999999999.9;
max[k] = -99999999999.9;
}


for(i=0; i<sys.n_atoms; i++)
for(k=0; k<3; k++)
{
if(sys.atoms[i].coor[k] < min[k]) min[k] = sys.atoms[i].coor[k];
if(sys.atoms[i].coor[k] > max[k]) max[k] = sys.atoms[i].coor[k];
}
//printf("max = %lf %lf %lf min =  %lf %lf %lf\n", 
//max[0], max[1], max[2],
//min[0], min[1], min[2]);

rave = 0.5 * sqrt( 
       (max[0] - min[0])*(max[0] - min[0]) 
     + (max[1] - min[1])*(max[1] - min[1]) 
     + (max[2] - min[2])*(max[2] - min[2]) 
       )/sqrt(3.0);
if(rave < 1.0) rave = 1.0;
//printf("r_ave = %lf\n", rave);
if(flag_par.msms)
	srf=&srf_ses;
else
	srf=&srf_sas;

for(i=0; i<sys.n_atoms; i++)
{
for(k=0; k <(*srf).at_srf[i].n_isrf; k++) // prima il contributo dalla propria(*srf)
{
l =(*srf).at_srf[i].isrf[k];
if(l >(*srf).n_srf_pt_tot) {printf("%i %i ....exiting....\n", l,(*srf).n_srf_pt_tot); exit(0);}
//printf("atom: %i surf nei %i pt. %i\n",i, k, l);
diffv(v,(*srf).srf_pt[l].r, sys.atoms[i].coor);
tmpf = modv(v)/rave; //distance scaled by rave 
if(tmpf > 0)
 {
 t0 =(*srf).srf_pt[l].a*dotv(v,(*srf).srf_pt[l].vec)/(tmpf*tmpf*tmpf);
 t1 = t0/tmpf;
 t2 = t1/tmpf;
 t3 = t2/tmpf;

 phi[i] = phi[i] + t0;
 ir1[i] = ir1[i] + t1; 
 ir2[i] = ir2[i] + t2; 
 ir3[i] = ir3[i] + t3; 
// if(i==742)
// printf("self: %lf %lf %lf %lf %lf\n", t0,t1,t2,t3, modv(v));
 }
// else {printf("Come fa ad essere <= 0???\n"); exit(0);}
}

for(j=0; j< neighbours[i].n_neighbours; j++) // poi dalle altre
for(k=0; k <(*srf).at_srf[neighbours[i].list[j]].n_isrf; k++)
{
l =(*srf).at_srf[neighbours[i].list[j]].isrf[k];
diffv(v,(*srf).srf_pt[l].r, sys.atoms[i].coor);
tmpf = modv(v)/rave; //distance scaled by rave 
 if(tmpf > 0)
 {
 t0 =(*srf).srf_pt[l].a*dotv(v,(*srf).srf_pt[l].vec)/(tmpf*tmpf*tmpf);
 t1 = t0/tmpf;
 t2 = t1/tmpf;
 t3 = t2/tmpf;
 phi[i] = phi[i] + t0;
 ir1[i] = ir1[i] + t1; 
 ir2[i] = ir2[i] + t2; 
 ir3[i] = ir3[i] + t3; 
// if(i==742)
// printf("oths: %lf %lf %lf %lf %lf\n", t0,t1,t2,t3, modv(v));
 }
// else {printf("Come fa ad essere <= 0???\n"); exit(0);}
}

//printf("atom %i: ir3[i] = %lf\n", i, ir3[i]);
}
for(i=0; i<sys.n_atoms; i++)
{
t0 = rave*rave*rave;
t1 = rave*t0;
t2 = rave*t1;
t3 = rave*t2;
if(phi[i] < 2*M_PI) phi[i] = 4*M_PI;
phi[i] = phi[i]/t0;
ir1[i] = ir1[i]/t1;
ir2[i] = ir2[i]/t2;
ir3[i] = ir3[i]/t3;
//phi[i] = 4 * M_PI;
if(ir1[i] > 0.0)
{
ir1[i] = phi[i]/ir1[i];
//printf("atom %i: ir = %lf %lf %lf\n", i, ir1[i], ir2[i], ir3[i]);
}
else {
    ir1[i] = phi[i]/abs(ir1[i]);
      printf("Unexpected negative radius at atom  %5i %4s %4s %4i %c\n atom radius set to average radius, it may be due to approximate numerical treatment of internal cavities or surfaces\nThis will be treated approximately\n", 
               sys.atoms[i].at_n,
               sys.atoms[i].at_name,
               sys.atoms[i].res_name,
               sys.atoms[i].res_n,
               sys.atoms[i].chain);
     }
if(ir2[i] > 0) 
ir2[i] = pow(phi[i]/ir2[i], 1.0/2.0);
//else ir2[i] = ir1[i] + dx;
if(ir3[i] > 0) 
ir3[i] = pow(phi[i]/ir3[i], 1.0/3.0);
//else ir3[i] = ir2[i] + dx;
if(ir2[i] > ir1[i] && ir1[i] > 0.0) ir2[i] = ir1[i];
if(ir3[i] > ir2[i] && ir2[i] > 0.0) ir3[i] = ir2[i];
if(flag_par.msms)
	gbr6 = gbr6_ses;
	else
	gbr6 = gbr6_sas;

//if(sys.n_atoms > RAVE_N_ATOMS)
//{
//if(ir3[i] < rave)
//gbr6[i] = ir3[i];
//else
//gbr6[i] = rave;
//}
//else gbr6[i] = ir3[i];
if(ir3[i] > 0) 
gbr6[i] = ir3[i];
else 
gbr6[i] = -1.0;
}

if(!flag_par.msms)
for(i=0; i<sys.n_atoms; i++)
{
if(gbr6[i] > 0.0)
{
if(gbr6[i] < sys.atoms[i].radius + flag_par.probe_radius) gbr6[i] = sys.atoms[i].radius +  flag_par.probe_radius;
if(gbr6[i] > flag_par.cutoff + flag_par.probe_radius) gbr6_ses[i] = flag_par.cutoff+ flag_par.probe_radius;
gbr6_ses[i] = (0.77322 + 0.06065 * flag_par.probe_radius) * gbr6[i] + (0.5142 -1.7585 * flag_par.probe_radius);
if(gbr6_ses[i] < sys.atoms[i].radius) gbr6_ses[i] = sys.atoms[i].radius;
if(gbr6_ses[i] > flag_par.cutoff ) gbr6_ses[i] = flag_par.cutoff;
}
else
{
gbr6[i] = DBL_MAX;
for(k=0; k< neighbours[i].n_neighbours; k++) // poi dalle altre
{
j = neighbours[i].list[k];
if(gbr6[j] > 0.0 && gbr6[j] < gbr6[i]) 
gbr6[i] = gbr6[j];
}
}
}
else
{
	srf = &srf_sas;
	gbr6 = gbr6_sas;
for(i=0; i<sys.n_atoms; i++)
{
 phi[i] = 0.0;
 ir1[i] = 0.0; 
 ir2[i] = 0.0; 
 ir3[i] = 0.0; 
for(k=0; k <(*srf).at_srf[i].n_isrf; k++) // prima il contributo dalla propria(*srf)
{
l =(*srf).at_srf[i].isrf[k];
if(l >(*srf).n_srf_pt_tot) {printf("%i %i ....exiting....\n", l,(*srf).n_srf_pt_tot); exit(0);}
//printf("atom: %i surf nei %i pt. %i\n",i, k, l);
diffv(v,(*srf).srf_pt[l].r, sys.atoms[i].coor);
//printf("%lf %lf %lf %lf\n", v[0],v[1],v[2], rave);
tmpf = modv(v)/rave; //distance scaled by rave 
 t0 =(*srf).srf_pt[l].a*dotv(v,(*srf).srf_pt[l].vec)/(tmpf*tmpf*tmpf);
 t1 = t0/tmpf;
 t2 = t1/tmpf;
 t3 = t2/tmpf;
//printf("%lf %lf %lf %lf\n", phi[i],ir1[i],ir2[i],ir3[i]);
//printf("%lf %lf %lf %lf\n", t0,t1,t2,t3);

 phi[i] = phi[i] + t0;
 ir1[i] = ir1[i] + t1; 
 ir2[i] = ir2[i] + t2; 
 ir3[i] = ir3[i] + t3; 
}
for(j=0; j< neighbours[i].n_neighbours; j++) // poi dalle altre
for(k=0; k <(*srf).at_srf[neighbours[i].list[j]].n_isrf; k++)
{
l =(*srf).at_srf[neighbours[i].list[j]].isrf[k];
diffv(v,(*srf).srf_pt[l].r, sys.atoms[i].coor);
tmpf = modv(v)/rave; //distance scaled by rave 
 t0 =(*srf).srf_pt[l].a*dotv(v,(*srf).srf_pt[l].vec)/(tmpf*tmpf*tmpf);
 t1 = t0/tmpf;
 t2 = t1/tmpf;
 t3 = t2/tmpf;
 phi[i] = phi[i] + t0;
 ir1[i] = ir1[i] + t1; 
 ir2[i] = ir2[i] + t2; 
 ir3[i] = ir3[i] + t3; 
}

//printf("atom %i: ir3[i] = %lf\n", i, ir3[i]);
}
for(i=0; i<sys.n_atoms; i++)
{
t0 = rave*rave*rave;
t1 = rave*t0;
t2 = rave*t1;
t3 = rave*t2;
phi[i] = phi[i]/t0;
ir1[i] = ir1[i]/t1;
ir2[i] = ir2[i]/t2;
ir3[i] = ir3[i]/t3;
//phi[i] = 4 * M_PI;
if(ir1[i] > 0.0)
{
ir1[i] = phi[i]/ir1[i];
//printf("atom %i: ir = %lf %lf %lf\n", i, ir1[i], ir2[i], ir3[i]);
}
else {
      ir1[i] = rave;
      printf("Unexpected negative radius at atom  %5i %4s %4s %4i %c \n atom radius set to average radius, it may be due to approximate numerical treatment of internal cavities\nTry with msms surface\n", 
               sys.atoms[i].at_n,
               sys.atoms[i].at_name,
               sys.atoms[i].res_name,
               sys.atoms[i].res_n,
               sys.atoms[i].chain);
     }
if(ir2[i] > 0) 
ir2[i] = pow(phi[i]/ir2[i], 1.0/2.0);
//else ir2[i] = ir1[i] + dx;
if(ir3[i] > 0) 
ir3[i] = pow(phi[i]/ir3[i], 1.0/3.0);
//else ir3[i] = ir2[i] + dx;
if(ir2[i] > ir1[i] && ir1[i] > 0.0) ir2[i] = ir1[i];
if(ir3[i] > ir2[i] && ir2[i] > 0.0) ir3[i] = ir2[i];
//if(sys.n_atoms > RAVE_N_ATOMS)
//{
//if(ir3[i] < rave)
//gbr6[i] = ir3[i];
//else
//gbr6[i] = rave;
//}
//else gbr6[i] = ir3[i]; 
gbr6[i] = ir3[i];
}
}
}

void srfpot(struct System *sys, struct Srf srf, double *gbr6, struct Neighbour *neighbours, struct Flag_par flag_par)
{
int i, j, k, l, jj;
double d, f, kd, k_el;
double      q = 1.602176487E-19;
double      N_av = 6.02214179E23;
double      e0 = 8.8541878176E-12;
double      kb = 1.3806504E-23;
double      RT = 1.3806504 * 6.02214179 * flag_par.temp;

k_el = N_av * 1E10 * q * q /(4 * M_PI * e0);
kd = 1E-10 * sqrt((1000 * 2 * flag_par.ions * N_av * q * q )/
                         (flag_par.sdie * e0 * kb * flag_par.temp));
for(i=0; i<(*sys).n_atoms; i++)
{
(*sys).atoms[i].temp = 0.0;
for(k=0; k < srf.at_srf[i].n_isrf; k++) // prima il contributo dalla propria srf
{
l = srf.at_srf[i].isrf[k];
d = distv(srf.srf_pt[l].r, (*sys).atoms[i].coor);
f = fgb(gbr6[i],flag_par.test_charge_radius,  flag_par.probe_radius,  flag_par.salt_radius,  d, flag_par.pdie,  flag_par.sdie,  kd, flag_par.kp);
(*sys).atoms[i].temp = (*sys).atoms[i].temp + (*sys).atoms[i].charge/(d *flag_par.pdie); 
(*sys).atoms[i].temp = (*sys).atoms[i].temp + f * (*sys).atoms[i].charge; 
for(j=0; j< neighbours[i].n_neighbours; j++) // poi dalle altre
{
jj = neighbours[i].list[j];
d = distv(srf.srf_pt[l].r, (*sys).atoms[jj].coor);
f = fgb(gbr6[jj],flag_par.test_charge_radius,  flag_par.probe_radius,  flag_par.salt_radius,  d, flag_par.pdie,  flag_par.sdie,  kd, flag_par.kp);
(*sys).atoms[i].temp = (*sys).atoms[i].temp + (*sys).atoms[jj].charge/(d *flag_par.pdie); 
(*sys).atoms[i].temp = (*sys).atoms[i].temp + f * (*sys).atoms[jj].charge; 
}
}
if(srf.at_srf[i].n_isrf != 0)
{
(*sys).atoms[i].temp = (k_el/1E3) * (*sys).atoms[i].temp / (double) srf.at_srf[i].n_isrf;
(*sys).atoms[i].occ = 1.0;
}
else 
{
(*sys).atoms[i].temp = 0.0;
(*sys).atoms[i].occ = 0.0;
}
}
}

void scalev(double t, double *v)
{
    v[0] = t*v[0];
    v[1] = t*v[1];
    v[2] = t*v[2];
}

void get_srf(struct System sys, struct Atom_grid *patom_grid, struct Srf *psrf_ses, struct Srf *psrf_sas, struct Flag_par flag_par)
{
      char buf[2048];//, fn[120], *tmps;
      FILE *fp;
      int i,j,k,l,n_faces, max_n_srf_pt;
      struct Srf tmp_srf, srf;
      struct Atom_grid atom_grid;
      double d1, d2, d3, sp, At, A, t, r_max;
      struct Msms_face *faces;
      
      atom_grid = *patom_grid;

        for(i=0; i< sys.n_atoms; i++)
        if((sys.atoms[i].radius + flag_par.probe_radius) > r_max)
            r_max = sys.atoms[i].radius + flag_par.probe_radius;
       (*psrf_ses).at_srf = (struct At_srf *) calloc( (size_t) sys.n_atoms, sizeof(struct At_srf));
        for(i=0; i< sys.n_atoms; i++)
        (*psrf_ses).at_srf[i].n_isrf=0;

        (*psrf_sas).at_srf = (struct At_srf *) calloc( (size_t) sys.n_atoms , sizeof(struct At_srf));
        for(i=0; i< sys.n_atoms; i++)
        (*psrf_sas).at_srf[i].n_isrf=0;

if(flag_par.msms) 
{
printf("Calling msms for molecular surface calculation...\n");
printf("Make sure that msms is in your PATH...\n");
printf("and that the minimum radius is set as to avoid one atom included in one another (e.g. use option -mr 1.0)\n");
sprintf(buf,"%s.xyzr",flag_par.file_msms);
fp = file_open(buf,"w");
for(i = 0; i< sys.n_atoms; i++)
fprintf(fp,"%11.3f %11.3f %11.3f %11.3f\n", 
sys.atoms[i].coor[0],
sys.atoms[i].coor[1],
sys.atoms[i].coor[2],
sys.atoms[i].radius);
fclose(fp);
sprintf(buf,"msms -if %s.xyzr -of %s -af %s -density %lf -probe_radius %lf -no_rest > %s.msms_log\n",flag_par.file_msms, flag_par.file_msms, flag_par.file_msms, 1.0/flag_par.msms_area, flag_par.probe_radius, flag_par.file_msms);
printf("running... %s\n", buf);
system(buf);
sprintf(buf,"%s.vert",flag_par.file_msms);
fp = file_open(buf,"r");
printf(".....I am reading surface points and normal vectors in %s.vert\n", flag_par.file_msms);
do 
{
fgets(buf,120,fp);
} 
while (buf[0] == '#');
sscanf(buf,"%i", &(tmp_srf.n_srf_pt_tot));
tmp_srf.srf_pt= (struct Srf_pt *) calloc((size_t) tmp_srf.n_srf_pt_tot, sizeof(struct Srf_pt));
        for(i=0; i < tmp_srf.n_srf_pt_tot ; i++)
        {
        tmp_srf.srf_pt[i].r=(double *) calloc((size_t) 3, sizeof(double));
        tmp_srf.srf_pt[i].vec=(double *) calloc((size_t) 3, sizeof(double));
        }

k=0;
while(fgets(buf, 120, fp) != NULL)
{
//printf("%s\n", buf);
sscanf(buf,"%lf %lf %lf %lf %lf %lf %i %i %i\n",
&(tmp_srf.srf_pt[k].r[0]),
&(tmp_srf.srf_pt[k].r[1]),
&(tmp_srf.srf_pt[k].r[2]),
&(tmp_srf.srf_pt[k].vec[0]),
&(tmp_srf.srf_pt[k].vec[1]),
&(tmp_srf.srf_pt[k].vec[2]),
&(tmp_srf.srf_pt[k].atn), // temporary reassign after reading all srf pts.
&j,
&(tmp_srf.srf_pt[k].type));
tmp_srf.srf_pt[k].at = j - 1;
tmp_srf.srf_pt[k].a= flag_par.msms_area;
//tmp_srf.at_srf[j-1].isrf[tmp_srf.at_srf[j-1].n_isrf++] = k;
k++;
if(k > tmp_srf.n_srf_pt_tot) {printf("more points in file than stated in header...\nexiting...\n"); exit(0);}
}
fclose(fp);
tmp_srf.n_srf_pt_tot = k;
printf("I have read %i surface points in %s.vert\n", tmp_srf.n_srf_pt_tot, flag_par.file_msms);
sprintf(buf,"%s.face",flag_par.file_msms);
fp = file_open(buf,"r");
printf(".....I am reading surface faces in %s.face\n", flag_par.file_msms);
do 
{
fgets(buf,120,fp);
} 
while (buf[0] == '#');
sscanf(buf,"%i", &n_faces);
faces = (struct Msms_face *) calloc((size_t) n_faces, sizeof(struct Msms_face));
k=0;
while(fgets(buf, 120, fp) != NULL)
{
sscanf(buf,"%i %i %i %i %*i\n",
&(faces[k].n1),
&(faces[k].n2),
&(faces[k].n3),
&(faces[k].type)
);
faces[k].n1 = faces[k].n1 - 1;
faces[k].n2 = faces[k].n2 - 1;
faces[k].n3 = faces[k].n3 - 1;
k++;
if(k > n_faces) {printf("more points in file than stated in header...\nexiting...\n"); exit(0);}
}
n_faces = k;
printf("I have read %i surface faces in %s.face\n", n_faces, flag_par.file_msms);
j = 0;
for(k = 0; k< n_faces; k++)
{
// 1 2 e 3 uguali
if((tmp_srf.srf_pt[faces[k].n1].at == tmp_srf.srf_pt[faces[k].n2].at)  &&
(tmp_srf.srf_pt[faces[k].n1].at == tmp_srf.srf_pt[faces[k].n3].at) )
{
l = tmp_srf.srf_pt[faces[k].n1].at;
(*psrf_ses).at_srf[l].n_isrf++;
j++;
}
else if(tmp_srf.srf_pt[faces[k].n1].at == tmp_srf.srf_pt[faces[k].n2].at)
{ 
l = tmp_srf.srf_pt[faces[k].n1].at;
(*psrf_ses).at_srf[l].n_isrf++;
j++;
l = tmp_srf.srf_pt[faces[k].n3].at;
(*psrf_ses).at_srf[l].n_isrf++;
j++;
}
else if(tmp_srf.srf_pt[faces[k].n1].at == tmp_srf.srf_pt[faces[k].n3].at)
{ 
l = tmp_srf.srf_pt[faces[k].n1].at;
(*psrf_ses).at_srf[l].n_isrf++;
j++;
l = tmp_srf.srf_pt[faces[k].n2].at;
(*psrf_ses).at_srf[l].n_isrf++;
j++;
}
else if(tmp_srf.srf_pt[faces[k].n2].at == tmp_srf.srf_pt[faces[k].n3].at)
{ 
l = tmp_srf.srf_pt[faces[k].n1].at;
(*psrf_ses).at_srf[l].n_isrf++;
j++;
l = tmp_srf.srf_pt[faces[k].n2].at;
(*psrf_ses).at_srf[l].n_isrf++;
j++;
}
else
{
l = tmp_srf.srf_pt[faces[k].n1].at;
(*psrf_ses).at_srf[l].n_isrf++;
j++;
l = tmp_srf.srf_pt[faces[k].n2].at;
(*psrf_ses).at_srf[l].n_isrf++;
j++;
l = tmp_srf.srf_pt[faces[k].n3].at;
(*psrf_ses).at_srf[l].n_isrf++;
j++;
}
}
(*psrf_ses).n_srf_pt_tot = j;
        (*psrf_ses).srf_pt= (struct Srf_pt*) calloc((size_t) (*psrf_ses).n_srf_pt_tot, sizeof(struct Srf_pt));
        for(i=0; i < (*psrf_ses).n_srf_pt_tot ; i++)
         {
         (*psrf_ses).srf_pt[i].r=(double *) calloc((size_t) 3, sizeof(double));
         (*psrf_ses).srf_pt[i].vec=(double *) calloc((size_t) 3, sizeof(double));
         }

//        for(i=0; i < (*tmp_srf).n_srf_pt_tot ; i++)
//	if((*tmp_srf).srf_pt[i].type == 3 || (*tmp_srf).srf_pt[i].type == 2)
j = 0;
for(k = 0, At = 0.0; k< n_faces; k++)
{
for(i = 0; i < 3 ; i++)
{
(*psrf_ses).srf_pt[j].r[i] = 
 (tmp_srf.srf_pt[faces[k].n1].r[i] + 
  tmp_srf.srf_pt[faces[k].n2].r[i] + 
  tmp_srf.srf_pt[faces[k].n3].r[i] )/3.0; 
(*psrf_ses).srf_pt[j].vec[i] = 
 (tmp_srf.srf_pt[faces[k].n1].vec[i] + 
  tmp_srf.srf_pt[faces[k].n2].vec[i] + 
  tmp_srf.srf_pt[faces[k].n3].vec[i] ); 
}
t = modv((*psrf_ses).srf_pt[j].vec);
for(i = 0; i < 3 ; i++)
  (*psrf_ses).srf_pt[j].vec[i]  = (*psrf_ses).srf_pt[j].vec[i] / t; 
d1 = distv(tmp_srf.srf_pt[faces[k].n1].r, tmp_srf.srf_pt[faces[k].n2].r);
d2 = distv(tmp_srf.srf_pt[faces[k].n2].r, tmp_srf.srf_pt[faces[k].n3].r);
d3 = distv(tmp_srf.srf_pt[faces[k].n3].r, tmp_srf.srf_pt[faces[k].n1].r);
sp = (d1 + d2 + d3)/2.0;
A = sqrt(sp*(sp-d1)*(sp-d2)*(sp-d3));
  (*psrf_ses).srf_pt[j].a  = A; 
At = At + A;

// All three face vertices belong to the same atom
(*psrf_ses).srf_pt[j].at = tmp_srf.srf_pt[faces[k].n1].at;
if(tmp_srf.srf_pt[faces[k].n1].at == tmp_srf.srf_pt[faces[k].n2].at  &&
tmp_srf.srf_pt[faces[k].n1].at == tmp_srf.srf_pt[faces[k].n3].at ) 
{
(*psrf_ses).srf_pt[j].at = tmp_srf.srf_pt[faces[k].n1].at;
(*psrf_ses).srf_pt[j].a = A;
(*psrf_ses).srf_pt[j].type = faces[k].type;
if(faces[k].type != 3 )
(*psrf_ses).srf_pt[j].type = -1;
(*psrf_ses).srf_pt[j].type = -2;
j++;
}
// First and second face vertices belong to the same atom and third one no
else if(
   tmp_srf.srf_pt[faces[k].n1].at == tmp_srf.srf_pt[faces[k].n2].at &&
   tmp_srf.srf_pt[faces[k].n2].at != tmp_srf.srf_pt[faces[k].n3].at) 
{
(*psrf_ses).srf_pt[j].at = tmp_srf.srf_pt[faces[k].n1].at;
(*psrf_ses).srf_pt[j].a = 2*A/3;
(*psrf_ses).srf_pt[j].type = faces[k].type;
if(faces[k].type != 3)
(*psrf_ses).srf_pt[j].type = -1;
j++;
(*psrf_ses).srf_pt[j].a = A/3;
(*psrf_ses).srf_pt[j].at = tmp_srf.srf_pt[faces[k].n3].at;
for(i = 0; i < 3 ; i++)
  {
  (*psrf_ses).srf_pt[j].r[i] = (*psrf_ses).srf_pt[j-1].r[i];
  (*psrf_ses).srf_pt[j].vec[i] = (*psrf_ses).srf_pt[j-1].vec[i];
  }
(*psrf_ses).srf_pt[j].type = faces[k].type ;
if(faces[k].type != 3)
(*psrf_ses).srf_pt[j].type = -1;
j++;
}
// First and third face vertices belong to the same atom and second one no
else if(
   tmp_srf.srf_pt[faces[k].n1].at != tmp_srf.srf_pt[faces[k].n2].at &&
   tmp_srf.srf_pt[faces[k].n1].at == tmp_srf.srf_pt[faces[k].n3].at) 
{
(*psrf_ses).srf_pt[j].a = 2*A/3;
(*psrf_ses).srf_pt[j].type = faces[k].type;
if(faces[k].type != 3)
(*psrf_ses).srf_pt[j].type = -1;
j++;
for(i = 0; i < 3 ; i++)
  {
  (*psrf_ses).srf_pt[j].r[i] = (*psrf_ses).srf_pt[j-1].r[i];
  (*psrf_ses).srf_pt[j].vec[i] = (*psrf_ses).srf_pt[j-1].vec[i];
  }
(*psrf_ses).srf_pt[j].a = A/3;
(*psrf_ses).srf_pt[j].at = tmp_srf.srf_pt[faces[k].n2].at;
(*psrf_ses).srf_pt[j].type = faces[k].type;
if(faces[k].type != 3)
(*psrf_ses).srf_pt[j].type = -1;
j++;
}
// Second and third face vertices belong to the same atom and first one no
else if(
   tmp_srf.srf_pt[faces[k].n1].at != tmp_srf.srf_pt[faces[k].n2].at &&
   tmp_srf.srf_pt[faces[k].n2].at == tmp_srf.srf_pt[faces[k].n3].at) 
{
(*psrf_ses).srf_pt[j].a = A/3;
(*psrf_ses).srf_pt[j].type = faces[k].type;
if(faces[k].type != 3)
(*psrf_ses).srf_pt[j].type = -1;
j++;
for(i = 0; i < 3 ; i++)
  {
  (*psrf_ses).srf_pt[j].r[i] = (*psrf_ses).srf_pt[j-1].r[i];
  (*psrf_ses).srf_pt[j].vec[i] = (*psrf_ses).srf_pt[j-1].vec[i];
  }
(*psrf_ses).srf_pt[j].a = 2*A/3;
(*psrf_ses).srf_pt[j].at = tmp_srf.srf_pt[faces[k].n2].at;
(*psrf_ses).srf_pt[j].type = faces[k].type;
if(faces[k].type != 3)
(*psrf_ses).srf_pt[j].type = -1;
j++;
}
// All three face vertices belong to different atoms 
else if(
   tmp_srf.srf_pt[faces[k].n1].at != tmp_srf.srf_pt[faces[k].n2].at &&
   tmp_srf.srf_pt[faces[k].n1].at != tmp_srf.srf_pt[faces[k].n3].at &&
   tmp_srf.srf_pt[faces[k].n2].at != tmp_srf.srf_pt[faces[k].n3].at) 
{
(*psrf_ses).srf_pt[j].a = A/3;
(*psrf_ses).srf_pt[j].type = faces[k].type;
if(faces[k].type != 3)
(*psrf_ses).srf_pt[j].type = -1;
j++;
(*psrf_ses).srf_pt[j].a = A/3;
(*psrf_ses).srf_pt[j].at = tmp_srf.srf_pt[faces[k].n2].at;
for(i = 0; i < 3 ; i++)
  {
  (*psrf_ses).srf_pt[j].r[i] = (*psrf_ses).srf_pt[j-1].r[i];
  (*psrf_ses).srf_pt[j].vec[i] = (*psrf_ses).srf_pt[j-1].vec[i];
  }
(*psrf_ses).srf_pt[j].type = faces[k].type;
if(faces[k].type != 3)
(*psrf_ses).srf_pt[j].type = -1;
j++;
(*psrf_ses).srf_pt[j].a = A/3;
(*psrf_ses).srf_pt[j].at = tmp_srf.srf_pt[faces[k].n3].at;
for(i = 0; i < 3 ; i++)
  {
  (*psrf_ses).srf_pt[j].r[i] = (*psrf_ses).srf_pt[j-1].r[i];
  (*psrf_ses).srf_pt[j].vec[i] = (*psrf_ses).srf_pt[j-1].vec[i];
  }
(*psrf_ses).srf_pt[j].type = faces[k].type;
if(faces[k].type != 3)
(*psrf_ses).srf_pt[j].type = -1;
j++;
}
}

for(i = 0; i <  sys.n_atoms ; i++)
{
(*psrf_ses).at_srf[i].isrf = (int *) calloc((size_t) ((*psrf_ses).at_srf[i].n_isrf), sizeof(int));
(*psrf_ses).at_srf[i].n_isrf = 0;
}

for(i = 0; i <  (*psrf_ses).n_srf_pt_tot ; i++)
 (*psrf_ses).at_srf[(*psrf_ses).srf_pt[i].at].isrf[(*psrf_ses).at_srf[(*psrf_ses).srf_pt[i].at].n_isrf++] = i;
        for(i=0, A=0; i < (*psrf_ses).n_srf_pt_tot ; i++)
        A = (*psrf_ses).srf_pt[i].a + A;
       printf("Total SES area: %lf A^2\n", A);
}
//else
{
//printf(".....I am computing surface points and normal vectors\n");
max_n_srf_pt = (int) (sys.n_atoms * 4.0 * M_PI * r_max * r_max/ flag_par.msms_area );
//printf("%i\n", max_n_srf_pt);

        (*psrf_sas).srf_pt= (struct Srf_pt *) calloc((size_t) max_n_srf_pt, sizeof(struct Srf_pt));

        for(i=0; i < max_n_srf_pt ; i++)
         {
         (*psrf_sas).srf_pt[i].r=(double *) calloc((size_t) 3, sizeof(double));
         (*psrf_sas).srf_pt[i].vec=(double *) calloc((size_t) 3, sizeof(double));
         }
        for(i=0; i< sys.n_atoms; i++)
        (*psrf_sas).at_srf[i].isrf= (int *) calloc( (size_t) (max_n_srf_pt/ sys.n_atoms) , sizeof(int));
        (*psrf_sas).alloc_srf_pt = max_n_srf_pt;
//atoms_on_grid(&atom_grid, sys);
//atom_grid_realloc(&atom_grid); 
grid_srf_dens(sys, atom_grid, flag_par, psrf_sas);
}
//        printf("--------------------------------\n");
}
void print_radii(struct System sys, struct Trj trj, double *gbr6_ses, double *gbr6_sas, struct Flag_par flag_par)
{
	double *tmpr;
	int i;
	char *buf;
	FILE *fp;
buf= (char *) calloc((size_t) 1024, sizeof(char));
if(flag_par.pqg)
{
tmpr = (double *) calloc((size_t) sys.n_atoms, sizeof(double));
for(i=0; i<sys.n_atoms; i++)
	tmpr[i] = sys.atoms[i].radius;
for(i=0; i<sys.n_atoms; i++)
{
//printf("%lf %lf\n", gbr6[i], tmpr[i] );
if(gbr6_ses[i] < sys.atoms[i].radius)
	gbr6_ses[i] = sys.atoms[i].radius;
sys.atoms[i].occ = sys.atoms[i].charge;
sys.atoms[i].temp = gbr6_ses[i];
}
for(i = 0; i < sys.n_atoms; i++)
tmpr[i] = sys.atoms[i].radius;
strcpy(buf,flag_par.file_out);
strcat(buf,".pqg");
fp = fopen(buf,"w");
printf(".....I am writing pqg file in\n%s\n",buf);
write_PDB_atoms(fp, sys.n_atoms, sys.atoms, trj);
fclose(fp);
}
for(i=0; i<sys.n_atoms; i++)
sys.atoms[i].radius = tmpr[i]; //rimetti il raggio originale
for(i=0; i<sys.n_atoms; i++)
{
if(flag_par.msms) 
{
if( gbr6_ses[i] < sys.atoms[i].radius)
gbr6_ses[i] = sys.atoms[i].radius;
}
//gbr6 are gbr6 for ses now
else if( gbr6_sas[i] < (sys.atoms[i].radius + flag_par.probe_radius))
gbr6_sas[i] = sys.atoms[i].radius + flag_par.probe_radius;

//else if( gbr6[i] < (sys.atoms[i].radius + flag_par.probe_radius))
//gbr6[i] = sys.atoms[i].radius + flag_par.probe_radius;
}
// OUTPUT GBR6
strcpy(buf,flag_par.file_out);
strcat(buf,".gbr");
fp = fopen(buf,"w");
        printf(".....I am writing GB radii in:\n%s\n",buf);
fprintf(fp, "       N  ATOM  RES RESN CH.   GBR6(SES)  GBR6(SAS)\n");
for(i=0; i<sys.n_atoms; i++)
   fprintf(fp, "GBR %5i %4s %4s %4i %c  %9.3f  %9.3f\n", 
               sys.atoms[i].at_n,
               sys.atoms[i].at_name,
               sys.atoms[i].res_name,
               sys.atoms[i].res_n,
               sys.atoms[i].chain,
               gbr6_ses[i],
               gbr6_sas[i]
               );
fclose(fp);
        printf("--------------------------------\n");
free(tmpr);
free(buf);
}




/******* Inizio calcolo del potenziale alla superficie *******/
void print_srf_pot(struct System sys, struct Trj trj, double *gbr6, struct Srf srf_ses, struct Srf srf_sas, struct Neighbour *neighbours, struct Area_atom *area_atom, struct Flag_par flag_par)
{
      char buf[1024];
      FILE *fp;
      int i,j,jj;
      struct System model;
      double tmpf;
      char chain_old, res_ins_old, res_name_old[5];
      int res_n_old;
      double area;
      struct Srf *srf;
srf = &srf_sas;
if(flag_par.srfpot)
{    
//Output srfatpot
strcpy(buf,flag_par.file_out);
strcat(buf,".srfatpot");
fp = fopen(buf,"w");
printf(".....I am writing atomic surface potentials in\n%s\n",buf);
write_PDB_atoms(fp, sys.n_atoms, sys.atoms, trj);
fclose(fp);
        printf("--------------------------------\n");
// HERE
//system_alloc(&model);
model.atoms = (struct Atom *) calloc((size_t) sys.n_atoms, sizeof(struct Atom));
for(i=0; i< sys.n_atoms; i++)
area_atom[i].area = 0.0; 
for(i=0; i<(*srf).n_srf_pt_tot; i++)
area_atom[(*srf).srf_pt[i].at].area = area_atom[(*srf).srf_pt[i].at].area + (*srf).srf_pt[i].a;
for(i=0; i< sys.n_atoms; i++)
copy_atom(&(model.atoms[i]),sys.atoms[i]);
model.n_atoms = sys.n_atoms;

if(flag_par.srfpotav)
{
for(i=0; i< model.n_atoms; i++)
{
tmpf=area_atom[i].area; 
model.atoms[i].temp = sys.atoms[i].temp * area_atom[i].area;
for(jj=0; jj<neighbours[i].n_neighbours; jj++)
{
j = neighbours[i].list[jj];
if(distv(model.atoms[i].coor,model.atoms[j].coor) <=  flag_par.cutoff2 && i != j)
{
tmpf = tmpf + area_atom[j].area;
model.atoms[i].temp = model.atoms[i].temp + area_atom[j].area * model.atoms[j].temp;
}
}
if(area_atom[i].area > 0.0)
model.atoms[i].temp = model.atoms[i].temp/tmpf;
else
model.atoms[i].temp = 0.0;
}
strcpy(buf,flag_par.file_out);
strcat(buf,".srfatpotav");
fp = fopen(buf,"w");
printf(".....I am writing averaged atomic surface potentials in\n%s\n",buf);
write_PDB_atoms(fp, model.n_atoms, model.atoms, trj);
fclose(fp);
}
        printf("--------------------------------\n");
free(model.atoms);
}
      /* Se richiesto stampa la superficie */
if(flag_par.surface)
{
	if(flag_par.msms)
srf = &srf_ses;
	else
srf = &srf_sas;

strcpy(buf,flag_par.file_out);
strcat(buf,".srf");
fp = fopen(buf,"w");

printf(".....I am writing %i surface points in\n%s\n",(*srf).n_srf_pt_tot,buf);

for(i=0; i<(*srf).n_srf_pt_tot; i++)
fprintf(fp,"ATOM  %5i  SRF SRF S   1    %8.3lf%8.3lf%8.3lf  1.00%6.2lf\n", 
i,
(*srf).srf_pt[i].r[0] , 
(*srf).srf_pt[i].r[1] , 
(*srf).srf_pt[i].r[2], 0.0
);
fclose(fp);
        printf("--------------------------------\n");
strcpy(buf,flag_par.file_out);
strcat(buf,".at_r_n");
fp = fopen(buf,"w");
printf(".....I am writing surface point atom n., coordinate and normal versor in\n%s\n",buf);
fprintf(fp,"ATNUM     X        Y        Z        VX       VY       VZ   \n");
for(i=0; i<(*srf).n_srf_pt_tot; i++)
fprintf(fp,"%5i %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", 
sys.atoms[(*srf).srf_pt[i].at].at_n,
(*srf).srf_pt[i].r[0] , 
(*srf).srf_pt[i].r[1] , 
(*srf).srf_pt[i].r[2] ,
(*srf).srf_pt[i].vec[0] , 
(*srf).srf_pt[i].vec[1] , 
(*srf).srf_pt[i].vec[2] 
);
fclose(fp);
        printf("--------------------------------\n");


strcpy(buf,flag_par.file_out);
strcat(buf,".area");
printf(".....I am writing surface area by sys, chain, residue, atom in\n%s\n",buf);
fp = fopen(buf,"w");
//printf("n_srf_pt_tot = %i\n",(srf).n_srf_pt_tot);
for(i=0; i<(*srf).n_srf_pt_tot; i++)
area_atom[(*srf).srf_pt[i].at].area = area_atom[(*srf).srf_pt[i].at].area + (*srf).srf_pt[i].a; 
for(i=0; i<sys.n_atoms; i++)
{
strcpy(area_atom[i].at_name,sys.atoms[i].at_name);
area_atom[i].alt_loc = sys.atoms[i].alt_loc;
strcpy(area_atom[i].res_name, sys.atoms[i].res_name);
area_atom[i].chain = sys.atoms[i].chain;
if(area_atom[i].chain == ' ') area_atom[i].chain = '_';
area_atom[i].at_n = i+1;
area_atom[i].res_n = sys.atoms[i].res_n;
area_atom[i].res_ins = sys.atoms[i].res_ins;
}

qsort(area_atom, sys.n_atoms, sizeof(struct Area_atom), &cmp_area_atom_by_string);
//area per il sistema
for(i=0, area = 0.0; i<sys.n_atoms; i++)
area = area + area_atom[i].area;
fprintf(fp, "System area = %.5e A^2 \n", area);

//area per catena
for(i=0,area = 0.0; i<sys.n_atoms; i++)
{
if(i == 0 || area_atom[i].chain == chain_old)
area = area + area_atom[i].area;
else 
{
fprintf(fp, "Chain %c area = %.3e A^2 \n", chain_old, area);
area = area_atom[i].area;
}
chain_old = area_atom[i].chain;
if(i == sys.n_atoms - 1)
fprintf(fp, "Chain %c area = %.3e A^2 \n", chain_old, area);
}

//area per residuo
for(i=0,area = 0.0; i<sys.n_atoms; i++)
{
if(i == 0 || (area_atom[i].chain == chain_old && 
area_atom[i].res_n == res_n_old && area_atom[i].res_ins == res_ins_old))
area = area + area_atom[i].area;
else 
{
fprintf(fp, "Residue %4s %5i %c chain %c area = %.3e A^2 \n", res_name_old, res_n_old, res_ins_old, chain_old, area);
area = area_atom[i].area;
}
chain_old = area_atom[i].chain;
res_n_old = area_atom[i].res_n;
res_ins_old = area_atom[i].res_ins;
strcpy(res_name_old, area_atom[i].res_name);
if(i == sys.n_atoms - 1)
fprintf(fp, "Residue %4s %5i %c chain %c area = %.3e A^2 \n", res_name_old, res_n_old, res_ins_old, chain_old, area);
}
qsort(area_atom, sys.n_atoms, sizeof(struct Area_atom), &cmp_area_atom_by_n);
for(i=0; i<sys.n_atoms; i++)
fprintf(fp, "Atom %4s %4s %5i %c chain %c area = %.3e A^2 \n", area_atom[i].at_name, area_atom[i].res_name, area_atom[i].res_n, area_atom[i].res_ins, area_atom[i].chain, area_atom[i].area);

fclose(fp);

        printf("--------------------------------\n");
}
}
/***************** Fine output superficie ****************************/

/***************** Inizio calcolo energia elettrostatica e forze ***********/
void calc_nrg(struct System sys, struct Energy *pnrg, double *gbr6, struct Srf srf, struct Neighbour *neighbours, struct Const constants, struct Flag_par flag_par)
{
int i, j, jj;
double t1, d, *v, f, *tmpgbr;
struct Energy nrg;
v=(double *) calloc((size_t) 3, sizeof(double));

if(flag_par.solv_e) 
{
(*pnrg).solv=0; 
(*pnrg).coul=0; 
(*pnrg).area=0; 
(*pnrg).born = 0;
for(i=0; i<(srf).n_srf_pt_tot; i++)
(*pnrg).area = (*pnrg).area + (srf).srf_pt[i].a;

(*pnrg).area = (*pnrg).area * flag_par.gamma;
if(!flag_par.msms)
{
	tmpgbr = (double *) calloc((size_t) sys.n_atoms, sizeof(double));
for(i=0; i<sys.n_atoms; i++)
{
tmpgbr[i] = gbr6[i];
gbr6[i] = (0.77322 + 0.06065 * flag_par.probe_radius) * gbr6[i] + (0.5142 -1.7585 * flag_par.probe_radius);
if(gbr6[i] < sys.atoms[i].radius)
        gbr6[i] = sys.atoms[i].radius;
}
}
for(i=0; i<sys.n_atoms; i++)
{
for(jj=0; jj<neighbours[i].n_neighbours; jj++)
if(neighbours[i].list[jj] > i)
{
j = neighbours[i].list[jj];
t1 = constants.k_el * sys.atoms[i].charge * sys.atoms[j].charge/1000;
diffv(v,sys.atoms[i].coor, sys.atoms[j].coor);
d = modv(v);
f = fgb(gbr6[i],gbr6[j],  flag_par.probe_radius,  flag_par.salt_radius,  d, flag_par.pdie,  flag_par.sdie,  constants.kd, flag_par.kp);
(*pnrg).solv_a[i] = (*pnrg).solv_a[i] + t1 * f/2.0;
(*pnrg).solv_a[j] = (*pnrg).solv_a[j] + t1 * f/2.0;
(*pnrg).solv = (*pnrg).solv + t1 * f;
if(d<=0) {printf("Distance between atoms %i and %i < 0, exiting....\n", i, j); exit(0);}
(*pnrg).coul = (*pnrg).coul + t1/(flag_par.pdie * d);
}
// aggiungo il self
d = 0;
t1 = constants.k_el * sys.atoms[i].charge * sys.atoms[i].charge/1000.0;
//printf("%lf\n", fgb(gbr6[i],gbr6[i],  flag_par.probe_radius,  flag_par.salt_radius,  d, flag_par.pdie,  flag_par.sdie,  kd, flag_par.kp));
f = fgb(gbr6[i],gbr6[i],  flag_par.probe_radius,  flag_par.salt_radius,  d, flag_par.pdie,  flag_par.sdie,  constants.kd, flag_par.kp);
(*pnrg).born = (*pnrg).born + t1 * f / 2.0;
(*pnrg).solv_a[i] = (*pnrg).solv_a[i] + t1 * f / 2.0 ;
(*pnrg).solv = (*pnrg).solv + t1 * f / 2.0;
}
}
if(!flag_par.msms)
for(i=0; i<sys.n_atoms; i++)
gbr6[i] = tmpgbr[i];
}

void print_nrg(struct Energy nrg, struct System sys, struct Flag_par flag_par)
{
	FILE *fp;
	int i;
	char buf[1024];
printf("ENERGIES\n");
printf("--------------------------------\n");
printf("Total Coulomb      energy: %14.4lf (kJ/mol)\n", nrg.coul);
printf("Total el_solvation energy: %14.4lf (kJ/mol)\n", nrg.solv);
printf(" --Born self       energy: %14.4lf (kJ/mol)\n", nrg.born);
printf(" --Coulomb solv    energy: %14.4lf (kJ/mol)\n", nrg.solv - nrg.born);
printf("Total sasa         energy: %14.4lf (kJ/mol)\n", nrg.area);
printf("Total              energy: %14.4lf (kJ/mol)\n", (nrg.coul + nrg.solv + nrg.area));
printf("--------------------------------\n");
 
strcpy(buf,flag_par.file_out);
strcat(buf,".solv_nrg");
fp = fopen(buf,"w");
printf(".....I am writing total and atomic solvation energies in\n%s\n",buf);
fprintf(fp, "Total Coulomb      energy: %14.4lf (kJ/mol)\n", nrg.coul);
fprintf(fp, "Total el_solvation energy: %14.4lf (kJ/mol)\n", nrg.solv);
fprintf(fp, " --Born self       energy: %14.4lf (kJ/mol)\n", nrg.born);
fprintf(fp, " --Coulomb solv    energy: %14.4lf (kJ/mol)\n", nrg.solv - nrg.born);
fprintf(fp, "Total sasa         energy: %14.4lf (kJ/mol)\n", nrg.area);
fprintf(fp, "Total              energy: %14.4lf (kJ/mol)\n", (nrg.coul + nrg.solv + nrg.area));
   fprintf(fp, "           N   ATOM  RES RESN CH. CHARGE SOLVATION_ENERGY (kJ/mol)\n"); 
for(i=0; i<sys.n_atoms; i++)
   fprintf(fp, "SOLV NRG %5i %4s %4s %4i %c  %7.3f %12.3f\n", 
               sys.atoms[i].at_n,
               sys.atoms[i].at_name,
               sys.atoms[i].res_name,
               sys.atoms[i].res_n,
               sys.atoms[i].chain,
               sys.atoms[i].charge,
               nrg.solv_a[i]
               );
fclose(fp);
        printf("--------------------------------\n");
}

/***************** Inizio calcolo grid elettrostatica ***********/

void out_grid(struct System sys, struct Atom_grid atom_grid, double *gbr6, struct Neighbour *neighbours, struct Const constants, struct Flag_par flag_par)
{
	double xmin,xmax,ymin,ymax,zmin,zmax,t1;
	double hx, hy, hz;
	double f, d, *v;
        int i, j, k, l, ii, jj, kk, m;
	int nx, ny, nz,il,jl,kl;
	double ***out;
	bool ***isin;
	char buf[1024];
	FILE *fp;
	nx = flag_par.nx;
	ny = flag_par.ny;
	nz = flag_par.nz;
v = (double *) calloc((size_t) 3, sizeof(double));
if(flag_par.grid)
{
if(nx > 0 && ny > 0 && nz > 0)
{
if(1.0/constants.kd > 10.0) t1 = 0.1; else t1 = constants.kd;
xmin = atom_grid.x_min - 3.0/t1;
xmax = atom_grid.x_max + 3.0/t1;
ymin = atom_grid.y_min - 3.0/t1;
ymax = atom_grid.y_max + 3.0/t1;
zmin = atom_grid.z_min - 3.0/t1;
zmax = atom_grid.z_max + 3.0/t1;
hx = (xmax - xmin)/(double) nx;
hy = (ymax - ymin)/(double) ny;
hz = (zmax - zmin)/(double) nz;
}
else
{
xmin = (atom_grid.x_min + atom_grid.x_max)/2.0 - (double) nx * flag_par.mesh/2.0;
ymin = (atom_grid.y_min + atom_grid.y_max)/2.0 - (double) ny * flag_par.mesh/2.0;
zmin = (atom_grid.z_min + atom_grid.z_max)/2.0 - (double) nz * flag_par.mesh/2.0;
xmin = (atom_grid.x_min + atom_grid.x_max)/2.0 + (double) nx * flag_par.mesh/2.0;
ymin = (atom_grid.y_min + atom_grid.y_max)/2.0 + (double) ny * flag_par.mesh/2.0;
zmin = (atom_grid.z_min + atom_grid.z_max)/2.0 - (double) nz * flag_par.mesh/2.0;
hx = hy = hz = flag_par.mesh;
}

out = (double ***) calloc((size_t) nx, sizeof(double **));
for(i = 0; i < nx; i++)
{
out[i] = (double **) calloc((size_t) ny, sizeof(double *));
for(j = 0; j < ny; j++)
out[i][j] = (double *) calloc((size_t) nz, sizeof(double));
}
isin = (bool ***) calloc((size_t) nx, sizeof(bool **));
for(i = 0; i < nx; i++)
{
isin[i] = (bool **) calloc((size_t) ny, sizeof(bool *));
for(j = 0; j < ny; j++)
isin[i][j] = (bool *) calloc((size_t) nz, sizeof(bool));
}

strcpy(buf,flag_par.file_out);
strcat(buf,".dx");
fp = fopen(buf,"w");
printf(".....I am writing the potential in DX format in\n%s\n",buf);
fprintf(fp,"# Comments\n\
object 1 class gridpositions counts %i %i %i\n\
origin %f %f %f\n\
delta %f 0.0 0.0\n\
delta 0.0 %f 0.0\n\
delta 0.0 0.0 %f\n\
object 2 class gridconnections counts %i %i %i\n\
object 3 class array type double rank 0 items %i data follows\n",
nx, ny, nz, xmin, ymin, zmin,hx,hy,hz,nx,ny,nz, nx*ny*nz );
if( (flag_par.cutoff > 3.0/t1 - hx || flag_par.cutoff > 3.0/t1 - hy || flag_par.cutoff > 3.0/t1 - hz) && (nx == -1 || ny == -1 || nz == -1)) 
{
printf("Reset cutoff for grid calculation to %lf...\n", 3.0/t1);
il = (int) (3.0/t1)/hx - 1;
jl = (int) (3.0/t1)/hy - 1;
kl = (int) (3.0/t1)/hz - 1;
}
else
{
il = (int) (flag_par.cutoff)/hx + 1;
jl = (int) (flag_par.cutoff)/hy + 1;
kl = (int) (flag_par.cutoff)/hz + 1;
}
for(m=0; m<sys.n_atoms; m++)
{
i = (int) ((sys.atoms[m].coor[0] - xmin) / hx);
j = (int) ((sys.atoms[m].coor[1] - ymin) / hy);
k = (int) ((sys.atoms[m].coor[2] - zmin) / hz);
for(ii = -il; ii <= il ; ii++)
for(jj = -jl; jj <= jl ; jj++)
for(kk = -kl; kk <= kl ; kk++)
if(!isin[i+ii][j+jj][k+kk])
{
v[0] = xmin + (double) (i+ii) * hx;
v[1] = ymin + (double) (j+jj) * hy;
v[2] = zmin + (double) (k+kk) * hz;
d = distv( sys.atoms[m].coor, v);
f = fgb(gbr6[i],flag_par.test_charge_radius,  flag_par.probe_radius,  flag_par.salt_radius,  d, flag_par.pdie,  flag_par.sdie,  constants.kd, flag_par.kp);
out[i+ii][j+jj][k+kk] = out[i+ii][j+jj][k+kk] + sys.atoms[m].charge/(d *flag_par.pdie); 
out[i+ii][j+jj][k+kk] = out[i+ii][j+jj][k+kk] + sys.atoms[m].charge*f; 
if (d < (sys.atoms[m].radius + flag_par.probe_radius)) 
{
out[i+ii][j+jj][k+kk]=0.0;
isin[i+ii][j+jj][k+kk]=1;
}
}
}

for(i = 0,l=1; i < nx; i++)
for(j = 0; j < ny; j++)
for(k = 0; k < nz; k++)
{
out[i][j][k] = constants.k_el * out[i][j][k] /1000;
fprintf(fp, "%f ", out[i][j][k]);
if(l%3==0) fprintf(fp,"\n");
l++;
}
if( (l+1)%3 != 2) fprintf(fp,"\n");
fprintf(fp,"attribute \"dep\" string \"positions\"\n\
object \"regular positions regular connections\" class field\n\
component \"positions\" value 1\n\
component \"connections\" value 2\n\
component \"data\" value 3\n");
fclose(fp);
        printf("--------------------------------\n");
}
/********** Fine calcolo griglia ****************/
}

void calc_print_pka(struct System sys, struct Trj trj, double *gbr6, struct Neighbour *neighbours, struct Const constants, struct Flag_par *flag_par)
{
      int i,j,k,l,m,mm,n,jj,kk,tmpi,n_pka=14,isave,max_n_at_res;
      double d, **cova, f, t1, t2, t3, r_acc;
      FILE *fp;
      char buf[256];
      struct pka_def *pkadefs;
      struct System model, sites;
      double *energy, *charge_tot, charge, *gbr6m, *tmpgbr6;
      struct Neighbour *nei_sites;
      struct Srf *psrf_ses, srf_ses;
      struct Srf *psrf_sas, srf_sas;
      char pka_atom[][4]={"CD" ,"CG","ND1","ND1","NE2","NE2","NE2","NE2","CZ" ,"NZ", "OH", "SG","N",  "C"};
      char pka_res[][4] ={"GLU","ASP","HIE","HSE","HID","HSD","HIS","HIP","ARG","LYS","TYR","CYS","FRS","LST"};
//      double pka_pka[] ={4.4, 4.0, 6.3, 6.3, 12.0, 10.4, 9.6, 8.3, 7.5, 3.2};
//      double pka_pka[] ={4.5, 3.9, 6.8, 6.8, 13.8, 10.4, 10.2, 6.2, 7.5, 3.2};
      double pka_pka[] ={4.61, 3.97, 6.63, 6.63, 6.63, 6.63, 6.63, 6.63, 13.8, 10.63, 10.61, 6.25, 7.5, 3.2};
      double pka_charge[] ={-1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 1.0, -1.0};
      double pka_charge_add[] ={1.0, 1.0, 0.0, 0.0, 0.0, 0.0,-1.0, 0.0, -1.0, -1.0, 0.0, 0.0,-1.0,  1.0};
      double *pka, *pka0, **g, *bg, *bgm, *gm, *site_charge, *pka_self_p, *pka_self_m, *pka_eim, *pka_eip, *pka_em, *pka_ep, *dpka, pH, *state; 
      struct Atom_grid atom_grid; 
      double *ph, pH_beg = -4.0, pH_end = 18.0, dpH = 0.5, **titration, *pka_scale;
      int *pka_site, n_sites, *pka_res_index, flag;
//      double w1=1.53, w2=0.64, w3=0.69;
      double w1=0.96, w2=0.34, w3=0.47;
//      double w1=1, w2=1, w3=1;
      int pka_state[] ={1, 1, 1, 1, 1, 1, 0, 0, 1, 1}, *pka_pair;
      int nsteps, presteps, npH;

isave = (*flag_par).msms;
(*flag_par).msms=0;

if((*flag_par).pka)
{
/********** Inizio calcolo pka *****************/
/* 0) alloca memoria
   1) trova i siti ionizzabili assegna carica e pka
   2) calcola energie di background e funzioni di Green
      aggiunge una energia grande all Green function fra siti sulla stessa 
      istidina
   3) lancia il MC
   4) analizza il MC
*/
m = 0;
for(i=0; i< sys.n_atoms; i++)
sys.atoms[i].temp = gbr6[i];
make_system(&sys, &trj);
for(i=0; i< sys.n_atoms; i++)
gbr6[i] = sys.atoms[i].temp ;

max_n_at_res = 0;
for(i=0; i< sys.n_residues; i++)
{
j=sys.residues[i].end - sys.residues[i].beg + 1;
if(j > max_n_at_res) max_n_at_res = j;
}

if((*flag_par).defpka)
{
	fp=fopen((*flag_par).pkadef,"r");
	i=0;
        while(fgets(buf,120,fp) != NULL) i++;
	rewind(fp);
	n_pka=i;
	i=0;
        pkadefs=(struct pka_def *) calloc((size_t) n_pka,sizeof(struct pka_def));	
        while(fgets(buf,120,fp) != NULL) 
	{
        sscanf(buf,"%s %s %lf %i %lf", 
			pkadefs[i].res,
			pkadefs[i].atom,
			&(pkadefs[i].pka),
			&(pkadefs[i].state),
			&(pkadefs[i].charge)
			);
			pkadefs[i].charge_add = -pkadefs[i].charge * (double) pkadefs[i].state; 
	i++;
	}
}
else // assign known pkas
{
pkadefs=(struct pka_def *) calloc((size_t) n_pka,sizeof(struct pka_def));	
for(i=0; i<n_pka; i++)
{
			strcpy(pkadefs[i].res,pka_res[i]);
			strcpy(pkadefs[i].atom,pka_atom[i]);
			pkadefs[i].pka = pka_pka[i];
			pkadefs[i].state = pka_state[i];
			pkadefs[i].charge = pka_charge[i];
			pkadefs[i].charge_add = pka_charge_add[i];
}
}
/*
for(i=0; i< n_pka; i++)
        printf("%s %s %lf %i %lf %lf\n", 
			pkadefs[i].res,
			pkadefs[i].atom,
			(pkadefs[i].pka),
			(pkadefs[i].state),
			(pkadefs[i].charge),
			(pkadefs[i].charge_add)
			);
exit(0);
*/
/* conta i siti ionizzabili */
for(k=0; k< sys.n_chains; k++)
for(i=sys.chains[k].beg; i<= sys.chains[k].end; i++)
{
for(j=0; j< n_pka; j++)
if((!strcmp(sys.atoms[i].res_name,pkadefs[j].res) && 
    !strcmp(sys.atoms[i].at_name,pkadefs[j].atom)))
{ 
// exclude ssbonded cysteines 
tmpi = 1;
if(!strcmp("CYS",sys.atoms[i].res_name) &&
   !strcmp("SG",sys.atoms[i].at_name))
for(l=0; l< sys.n_residues;l++)
if(!strcmp(sys.residues[l].res_name,"CYS"))
for(n=sys.residues[l].beg; n<=sys.residues[l].end; n++)
if(!strcmp(sys.atoms[n].at_name,"SG") || !strcmp(sys.atoms[n].at_name,"SG1"))
{
d = distv(sys.atoms[i].coor,sys.atoms[n].coor);
if(d <= 3.0 && i != n) 
{
tmpi = 0;
strcpy(sys.atoms[i].at_name,"SX");
strcpy(sys.atoms[n].at_name,"SX");
}
}
if(tmpi)
m++;
}
if(!(*flag_par).defpka && (*flag_par).autoter)
{
   if((!strcmp(sys.atoms[i].at_name,"N") && 
      sys.atoms[i].res_n == sys.atoms[sys.chains[k].beg].res_n
      && sys.atoms[i].res_ins == sys.atoms[sys.chains[k].beg].res_ins)
   || (!strcmp(sys.atoms[i].at_name,"C") && 
      sys.atoms[i].res_n == sys.atoms[sys.chains[k].end].res_n
      && sys.atoms[i].res_ins == sys.atoms[sys.chains[k].end].res_ins)
     )
m++;
}
}
n_sites = m;
/* qui alloca */
pka0 = (double *) calloc((size_t) m, sizeof(double));
pka = (double *) calloc((size_t) m, sizeof(double));
bg = (double *) calloc((size_t) m, sizeof(double));
site_charge = (double *) calloc((size_t) m, sizeof(double));
pka_site = (int *) calloc((size_t) m, sizeof(int));
pka_res_index = (int *) calloc((size_t) m, sizeof(int));
state = (double *) calloc((size_t) m, sizeof(double));
g = (double **) calloc((size_t) m, sizeof(double *));
for(i=0;i<m; i++)
g[i] = (double *) calloc((size_t) m, sizeof(double));
cova = (double **) calloc((size_t) m, sizeof(double *));
for(i=0;i<m; i++)
cova[i] = (double *) calloc((size_t) m, sizeof(double));

gm = (double *) calloc((size_t) m, sizeof(double));
bgm = (double *) calloc((size_t) m, sizeof(double));
pka_self_p = (double *) calloc((size_t) m, sizeof(double)); // self 
pka_self_m = (double *) calloc((size_t) m, sizeof(double));
pka_eim = (double *) calloc((size_t) m, sizeof(double)); // bg
pka_eip = (double *) calloc((size_t) m, sizeof(double));
pka_em = (double *) calloc((size_t) m, sizeof(double)); // all
pka_ep = (double *) calloc((size_t) m, sizeof(double));
dpka = (double *) calloc((size_t) m, sizeof(double));
pka_pair = (int *) calloc((size_t) m, sizeof(int));
pka_scale = (double *) calloc((size_t) m, sizeof(double));
for(i=0; i< m; i++)
pka_pair[i] = -1;
titration = (double **) calloc((size_t) m, sizeof(double *));
npH = (int) ((pH_end - pH_beg)/dpH) + 1;
ph = (double *) calloc((size_t) npH, sizeof(double));
for(i=0; i< npH; i++)
ph[i] = pH_beg + (double) i * dpH;

for(i=0;i<m; i++)
titration[i] = (double *) calloc((size_t) npH, sizeof(double));
energy = (double *) calloc((size_t) npH, sizeof(double));
charge_tot = (double *) calloc((size_t) npH, sizeof(double));

printf("%i ionizable sites found\n", n_sites);

// assegna i valori imperturbati 
l = 0;
for(k=0; k< sys.n_chains; k++)
for(i=sys.chains[k].beg; i<= sys.chains[k].end; i++)
{
flag = 0;
for(j=0; j< n_pka; j++)
if((!strcmp(sys.atoms[i].res_name,pkadefs[j].res) && 
    !strcmp(sys.atoms[i].at_name,pkadefs[j].atom)))
 {
 site_charge[l] = pkadefs[j].charge;
 pka0[l] = pkadefs[j].pka;
 pka[l] = pkadefs[j].pka;
 pka_site[l] = i;
 sys.atoms[i].charge = sys.atoms[i].charge + pkadefs[j].charge_add;
 flag = 1;
 }
if(!(*flag_par).defpka && (*flag_par).autoter)
{
   if((!strcmp(sys.atoms[i].at_name,"N") && 
      sys.atoms[i].res_n == sys.atoms[sys.chains[k].beg].res_n
      && sys.atoms[i].res_ins == sys.atoms[sys.chains[k].beg].res_ins))
 {
 site_charge[l] = pkadefs[8].charge;
 pka0[l] = pkadefs[8].pka;
 pka[l] = pkadefs[8].pka;
 pka_site[l] = i;
 sys.atoms[i].charge = sys.atoms[i].charge + pkadefs[8].charge_add;
 flag = 1;
 }
   if((!strcmp(sys.atoms[i].at_name,"C") && 
      sys.atoms[i].res_n == sys.atoms[sys.chains[k].end].res_n
      && sys.atoms[i].res_ins == sys.atoms[sys.chains[k].end].res_ins)
     )
 {
 site_charge[l] = pkadefs[9].charge;
 pka0[l] = pkadefs[9].pka;
 pka[l] = pkadefs[9].pka;
 pka_site[l] = i;
 sys.atoms[i].charge = sys.atoms[i].charge + pkadefs[9].charge_add;
 flag = 1;
 }
}
 if(flag)
 {
 for(j=0; j < sys.n_residues; j++)
 if(i >= sys.residues[j].beg && i <= sys.residues[j].end)
 pka_res_index[l] = j;
 l++;
 }
}
charge = 0.0;
for(i=0; i< sys.n_atoms; i++)
charge = charge + sys.atoms[i].charge;
printf("Total charge after charge add %lf\n", charge);
sites.atoms = (struct Atom *) calloc((size_t) n_sites, sizeof(struct Atom));
for(l=0; l< n_sites; l++)
copy_atom(&(sites.atoms[l]),sys.atoms[pka_site[l]]);
sites.n_atoms = n_sites;
get_grid_parameters(sites, (*flag_par).probe_radius, (*flag_par).cutoff, &atom_grid, 'v');
//printf("%i %i %i %lf\n", atom_grid.grid_X, atom_grid.grid_Y, atom_grid.grid_Z, atom_grid.mesh);
atom_grid.PBC_ON = 0;
atom_grid_alloc(&atom_grid);
atoms_on_grid(&atom_grid, sites);
atom_grid_realloc(&atom_grid); 
nei_sites = neighbour_alloc(n_sites, atom_grid.max_atom_node);
if(nei_sites == NULL) 
{
printf("Could not allocate neighbours for titratable sites\nexiting...\n");
exit(0);
}
printf("n_sites: %i\n", n_sites);
grid_neighbour(atom_grid, sites, (*flag_par).cutoff, nei_sites,'v');
neighbour_realloc(nei_sites, n_sites); 
//for(l=0; l< n_sites; l++)
//printf("site %i: n_neighbours: %i\n", l, nei_sites[l].n_neighbours);
//inizio calcolo energie sulla proteina 
for(l=0; l< n_sites; l++)
{
i = pka_site[l];
pka_eip[l] = pka_ep[l] = 0.0;
d = 0;
f = fgb(gbr6[i],gbr6[i],  (*flag_par).probe_radius,  (*flag_par).salt_radius,  d, (*flag_par).pkadie,  (*flag_par).sdie,  constants.kd, (*flag_par).kp);
pka_self_p[l] = (constants.k_el * 1.0 * 1.0) * f/2.0;
//printf("%lf %lf\n", pka_self_p[l],f);
pka_eip[l] = constants.k_el * site_charge[l] * sys.atoms[i].charge * f;
// qui la lista dei neighbours non include se' stesso 
for(jj=0; jj<neighbours[i].n_neighbours; jj++)
 {
 j = neighbours[i].list[jj];
 d = distv( sys.atoms[i].coor, sys.atoms[j].coor);
 f = fgb(gbr6[i],gbr6[j],  (*flag_par).probe_radius,  (*flag_par).salt_radius,  d, (*flag_par).pkadie,  (*flag_par).sdie,  constants.kd, (*flag_par).kp);
if(i != j)
 {
pka_eip[l] = pka_eip[l] + constants.k_el * site_charge[l] * sys.atoms[j].charge * f;  
pka_eip[l] = pka_eip[l] + constants.k_el * site_charge[l] * sys.atoms[j].charge / ((*flag_par).pkadie * d);  
//printf("pka_eip f: %12.5lf with charge %12.5lf at x = %8.3lf\n", f, sys.atoms[j].charge, sys.atoms[j].coor[0]);
 }
 }
// Qui si puo' intervenire per pesare
/*
if(gbr6[i] <= 4.0) pka_scale[l] = 1.0; 
else if(gbr6[i] >= 7.0) pka_scale[l] = 0.0; 
else pka_scale[l] = (7.0 - gbr6[i])/3.0;
//pka_scale[l]=1.0;
//printf("%i (gbr6 = %lf): prima %lf ", l, gbr6[i], pka_self_p[l]); 
pka_self_p[l] = pka_self_p[l] * (2.0 + (1/.4 - 2.0) * pka_scale[l]);
//printf("dopo %lf\n", pka_self_p[l]); 
if(pka_eip[l] > 0.0)
t1 = exp(-pka_eip[l]/(2.0*RT));
else 
t1 = exp(pka_eip[l]/(2.0*RT));
//t1 = 1.0;
pka_eip[l] = pka_eip[l] * (1.0 + (t1 - 1.0) * pka_scale[l]);
pka_ep[l] = pka_self_p[l] + pka_eip[l] * t1;
*/
// tolgo per il debug lo scaling
pka_ep[l] = pka_self_p[l] + pka_eip[l];
}

// Qui si calcolano gli elementi di interazione fra sito l ed m senza il fattore 1/2 
for(l=0; l< n_sites; l++)
for(mm=0; mm< nei_sites[l].n_neighbours; mm++)
{
m = nei_sites[l].list[mm];
i = pka_site[l];
j = pka_site[m];
d = distv( sys.atoms[i].coor, sys.atoms[j].coor);
f = fgb(gbr6[i],gbr6[j],  (*flag_par).probe_radius,  (*flag_par).salt_radius,  d, (*flag_par).pkadie,  (*flag_par).sdie,  constants.kd, (*flag_par).kp);
if(l!=m)
{
t1 = constants.k_el * site_charge[l] * site_charge[m] * (f + 1.0/(d * (*flag_par).pkadie));
//* (2.0 + (1/4.0 - 2.0) * sqrt(pka_scale[l] * pka_scale[m]) );
g[l][m] = t1;
}
else
g[l][m] = 0.0;
if(
!strcmp(sys.atoms[i].res_name,"HIS") && !strcmp(sys.atoms[j].res_name,"HIS") && pka_res_index[l] == pka_res_index[m])
if(
(!strcmp(sys.atoms[i].at_name,"ND1") && !strcmp(sys.atoms[j].at_name,"NE2")) ||
(!strcmp(sys.atoms[i].at_name,"NE2") && !strcmp(sys.atoms[j].at_name,"ND1")) 
)
{
g[l][m] = g[l][m] + LARGE_ENERGY;
pka_pair[l] = m;
pka_pair[m] = l;
}
}
atom_grid_free(&atom_grid);
// fine calcolo energie sulla proteina 
// calcolo sui model compound 
neighbour_free(neighbours, sys.n_atoms);
gbr6m = (double *) calloc((size_t) max_n_at_res, sizeof(double));
tmpgbr6 = (double *) calloc((size_t) max_n_at_res, sizeof(double));
for(l=0; l< n_sites; l++)
{
i = pka_site[l];
j = pka_res_index[l];
model.n_atoms = sys.residues[j].end - sys.residues[j].beg + 1;
model.atoms = (struct Atom *) calloc((size_t) model.n_atoms, sizeof(struct Atom));
m = 0;
for(k=sys.residues[j].beg; k<= sys.residues[j].end; k++)
copy_atom(&model.atoms[m++], sys.atoms[k]);
model.n_atoms = m;
//system_realloc(&model);
get_grid_parameters(model, (*flag_par).probe_radius, (*flag_par).cutoff, &atom_grid, 'v');
atom_grid_alloc(&atom_grid);
atoms_on_grid(&atom_grid, model);
atom_grid_realloc(&atom_grid); 

//allocate memory for neighbour lists
neighbours = neighbour_alloc(model.n_atoms, atom_grid.max_atom_node);
grid_neighbour(atom_grid, model, (*flag_par).cutoff, neighbours,'v');
neighbour_realloc(neighbours, model.n_atoms); 
get_srf(model, &atom_grid, &srf_ses, &srf_sas, (*flag_par));
atom_grid_free(&atom_grid);
pqr2gbr6(model, srf_ses, srf_sas, (*flag_par), neighbours, tmpgbr6, gbr6m);
srf_free(&srf_sas,model.n_atoms); // srf_ses is not allocated here
neighbour_free(neighbours, model.n_atoms);
i = pka_site[l] - sys.residues[pka_res_index[l]].beg;
// printf("%s %s\n", model.atoms[i].at_name, model.atoms[i].res_name); 

pka_eim[l] = pka_em[l] = 0.0;
d = 0;
f = fgb(gbr6m[i],gbr6m[i],  (*flag_par).probe_radius,  (*flag_par).salt_radius,  d, (*flag_par).pkadie,  (*flag_par).sdie,  constants.kd, (*flag_par).kp);
pka_self_m[l] = (constants.k_el * 1.0 * 1.0) * f/2.0;
//printf("%lf %lf\n", pka_self_m[l],f);
pka_eim[l] = pka_eim[l] + constants.k_el * site_charge[l] * model.atoms[i].charge * f;
//if(l == 25)
//printf("i = %i j = %i pka_eim = %lf\n", i, i, pka_eim[l]);
for(j=0; j<model.n_atoms; j++)
if(i != j)
{
//printf("%lf %lf %lf\n", model.atoms[i].coor[0], model.atoms[i].coor[1], model.atoms[i].coor[2]);
d = distv( model.atoms[i].coor, model.atoms[j].coor);
 {
f = fgb(gbr6[i],gbr6[j],  (*flag_par).probe_radius,  (*flag_par).salt_radius,  d, (*flag_par).pkadie,  (*flag_par).sdie,  constants.kd, (*flag_par).kp);
pka_eim[l] = pka_eim[l] + constants.k_el * site_charge[l] * model.atoms[j].charge * f;
pka_eim[l] = pka_eim[l] + constants.k_el * site_charge[l] * model.atoms[j].charge / ((*flag_par).pkadie * d);
//printf("pka_eim f: %12.5lf with charge %12.5lf at x = %8.3lf\n", f, model.atoms[j].charge, model.atoms[j].coor[0]);
 }
}
// Qui si puo' intervenire per pesare
pka_em[l] = pka_self_m[l] + pka_eim[l];
//system_free(&model);
free(model.atoms);
}
/* printf("ATOM RESI  RESN ---    pka_m  dpka_p\n"); */
for(l=0; l<n_sites; l++)
{
/* qui il site_charge da solo il segno dello shift */
i = pka_site[l];
dpka[l] = -site_charge[l] * (pka_ep[l] - pka_em[l])/(constants.N_av * constants.kb * (*flag_par).temp * log(10));
//printf("%4s %4s %5i --- %8e %8e\n",
//sys.atoms[i].at_name,
//sys.atoms[i].res_name,
//sys.atoms[i].res_n,
//pka[l],
//dpka[l]);
}
for(l=0; l<n_sites; l++)
state[l] = 0;
presteps = n_sites * (*flag_par).pkass/10;
nsteps = n_sites * (*flag_par).pkass;
srand(10);
for(pH = pH_beg, k = 0; pH <= pH_end; pH=pH + dpH, k++) 
{
printf(".....Doing MC at pH = %e\n", pH);
mc(state, site_charge, pka, dpka, g, cova, n_sites, nei_sites, (*flag_par).temp, pH, nsteps, presteps, &(energy[k]), &r_acc);
charge_tot[k] = 0.0;
for(l=0; l<n_sites; l++)
{
//printf("%8.3lf\n", state[l]);
titration[l][k] = state[l];
charge_tot[k] = charge_tot[k] + state[l] * site_charge[l];
}
}
npH = k;

strcpy(buf,(*flag_par).file_out);
strcat(buf,".ddg");
fp = fopen(buf,"w");
        printf(".....I am writing the pH-dependent free energy of folding in:\n%s\n",buf);
fprintf(fp,"#   pH    DDG (kJ/mol) charge\n");
for(k = 0; k < npH ; k++) 
fprintf(fp,"%6.2lf    %7.2lf    %7.2lf\n", pH_beg + k * dpH, energy[k], charge_tot[k]);
fclose(fp);
        printf("--------------------------------\n");

strcpy(buf,(*flag_par).file_out);
strcat(buf,".titration");
fp = fopen(buf,"w");
printf(".....I am writing the ionization state of all atoms at each pH in:\n%s\n",buf);
fprintf(fp, "ATOM  RES RESN CH.     pH  ionization\n"); 
for(l=0; l<n_sites; l++)
for(k=0; k < npH; k++)
{
i =  pka_site[l];
fprintf(fp,"%4s %4s %5i %c %8.3lf %8.3lf\n",sys.atoms[i].at_name, sys.atoms[i].res_name, sys.atoms[i].res_n, sys.atoms[i].chain, pH_beg + k * dpH, titration[l][k]);
}
fclose(fp);
        printf("--------------------------------\n");

for(l=0; l<n_sites; l++)
if(site_charge[l] < 0.0)  pka[l] = pH_beg - 1.0;
else if(site_charge[l] > 0.0)  pka[l] = pH_end + 1.0;


for(l=0; l<n_sites; l++)
{
if(pka_pair[l] == -1) 
for(k=0; k < npH-1; k++)
{
if(
   (titration[l][k] <= 0.5 && titration[l][k+1] > 0.5) ||
   (titration[l][k] >= 0.5 && titration[l][k+1] < 0.5) 
  )
pka[l] = ((0.5 - titration[l][k])/(titration[l][k+1] - titration[l][k])) * dpH + pH_beg + (double) k * dpH;
}
else
for(k=0; k < npH-1; k++)
{
m = pka_pair[l];
t1 = titration[l][k] + titration[m][k]; 
t2 = titration[l][k+1] + titration[m][k+1]; 
if(
    (t1 <= 0.5 &&  t2 > 0.5) ||
    (t1 >= 0.5 && t2 < 0.5) 
  )
pka[l] = ((0.5 - t1)/(t2 - t1)) * dpH + pH_beg + (double) k * dpH;
}
}

strcpy(buf,(*flag_par).file_out);
strcat(buf,".pka");
fp = fopen(buf,"w");
printf(".....I am writing pKas and dpKa contributions in:\n%s\n",buf);
fprintf(fp, "ATOM  RES RESN CH.     pKa    pKa_0 dpKa^self dpKa^bg  dpKa^ii    GBR6\n"); 
for(l=0; l<n_sites; l++)
{
t1 = -site_charge[l] * (pka_self_p[l] - pka_self_m[l])/(constants.N_av * constants.kb * (*flag_par).temp * log(10)) ; 
t2 = -site_charge[l] * (pka_eip[l] - pka_eim[l])/(constants.N_av * constants.kb * (*flag_par).temp * log(10)) ; 
t3 = pka[l] -pka0[l] -(t1 + t2);
//t2 = (dpka[l] - t1);
//t2 = t2 * 5.0/(5.0 + fabs(t2));
//t3 = (pka[l] - pka0[l] - dpka[l]);
//t3 = t3 * 3.0/(3.0 + fabs(t3));

//t1 = t1 * 3.0/(3.0 + fabs(t1)); 
//t2 = -site_charge[l] * (pka_eip[l] - pka_eim[l])/(N_av * kb * (*flag_par).temp * log(10)) ; 
//printf("%i %lf %lf %lf -- %lf\n", l, pka_eip[l], pka_eim[l], (pka_eip[l] - pka_eim[l]), -site_charge[l] * (pka_eip[l] - pka_eim[l])/(N_av * kb * flag_par.temp * log(10)));
/*
if(t2 > 1.0) t2 = 1.0;
if(t2 < -1.0) t2 = -1.0;
if(t3 > 1.0) t3 = 1.0;
if(t3 < -1.0) t3 = -1.0;
*/
i = pka_site[l];
t1 = t1 * 6.0/(6.0 + fabs(t1));
t2 = t2 * 6.0/(6.0 + fabs(t2));
t3 = t3 * 6.0/(6.0 + fabs(t3));
if(pka_pair[l] != -1)
{
if(l < pka_pair[l])
if(titration[l][0] < titration[pka_pair[l]][0])
 {
fprintf(fp,"%4s %4s %5i %c %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf\n",sys.atoms[i].at_name, sys.atoms[i].res_name, sys.atoms[i].res_n, sys.atoms[i].chain, 
//pka[l],
w1 * t1 + w2 * t2 + w3 * t3 + pka0[l],
pka0[l], t1, t2, t3,
gbr6[i]);
//pka_scale[l]);
 }
else
{
i = pka_site[pka_pair[l]];
fprintf(fp,"%4s %4s %5i %c %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf\n",sys.atoms[i].at_name, sys.atoms[i].res_name, sys.atoms[i].res_n, sys.atoms[i].chain, 
w1 * t1 + w2 * t2 + w3 * t3 + pka0[l],
pka0[l], t1, t2, t3, gbr6[i]);
}
}
else
 {
fprintf(fp,"%4s %4s %5i %c %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf\n",sys.atoms[i].at_name, sys.atoms[i].res_name, sys.atoms[i].res_n, sys.atoms[i].chain, 
w1 * t1 + w2 * t2 + w3 * t3 + pka0[l],
pka0[l], t1, t2, t3,gbr6[i]);
 }
}
fclose(fp);
        printf("--------------------------------\n");
}
/********** Fine calcolo pka *****************/
(*flag_par).msms=isave;
}


void get_grid_parameters(struct System sistema, double probe_radius, 
		double cutoff, struct Atom_grid *atom_grid, char type)
{
/* la logica e' di avere una griglia in cui due particelle allocate
su nodi non adiacenti siano sicuramente piu' lontani di 
se type e' 'v'
vdw1 + vdw2 + cutoff
vdw1 + vdw2 +  2 * probe_radius 
se type e' 'd' di
cutoff
*/
	int i;
	double safe = 0.1, max_vdw_ctc = 0.0, r_min, min_vol_atom;
        for(i=0; i< sistema.n_atoms; i++)
        if((2.0 * sistema.atoms[i].radius) > max_vdw_ctc)
            max_vdw_ctc = 2.0 * sistema.atoms[i].radius;
        r_min = max_vdw_ctc/2.0; 
        for(i=0; i< sistema.n_atoms; i++)
         if((sistema.atoms[i].radius) < r_min)
          r_min = sistema.atoms[i].radius;
        min_vol_atom = 4 * M_PI * r_min * r_min;

if(type='v')
{
if (probe_radius >= (cutoff/2.0) ) 
	(*atom_grid).mesh = max_vdw_ctc + 2*probe_radius + safe;
else
	(*atom_grid).mesh = max_vdw_ctc + cutoff + safe;
}
//printf("%lf %lf %lf mesh = %lf\n", max_vdw_ctc, probe_radius, cutoff, (*atom_grid).mesh);
if(type=='d')
(*atom_grid).mesh = cutoff + safe;

for (i=0; i< sistema.n_atoms; i++)
{
	 if(i == 0)
	     {
	       (*atom_grid).x_min = (*atom_grid).x_max = sistema.atoms[i].coor[0];
	       (*atom_grid).y_min = (*atom_grid).y_max = sistema.atoms[i].coor[1];
               (*atom_grid).z_min = (*atom_grid).z_max = sistema.atoms[i].coor[2];
              }
               if(sistema.atoms[i].coor[0] > (*atom_grid).x_max) (*atom_grid).x_max = sistema.atoms[i].coor[0];
               if(sistema.atoms[i].coor[1] > (*atom_grid).y_max) (*atom_grid).y_max = sistema.atoms[i].coor[1];
               if(sistema.atoms[i].coor[2] > (*atom_grid).z_max) (*atom_grid).z_max = sistema.atoms[i].coor[2];
               if(sistema.atoms[i].coor[0] < (*atom_grid).x_min) (*atom_grid).x_min = sistema.atoms[i].coor[0];
               if(sistema.atoms[i].coor[1] < (*atom_grid).y_min) (*atom_grid).y_min = sistema.atoms[i].coor[1];
               if(sistema.atoms[i].coor[2] < (*atom_grid).z_min) (*atom_grid).z_min = sistema.atoms[i].coor[2];
//	       printf("%lf %lf %lf\n", (*atom_grid).x_min, sistema.atoms[i].coor[0], (*atom_grid).x_max);
}
(*atom_grid).x_min = (*atom_grid).x_min - safe;
(*atom_grid).y_min = (*atom_grid).y_min - safe;
(*atom_grid).z_min = (*atom_grid).z_min - safe;
(*atom_grid).x_max = (*atom_grid).x_max + safe;
(*atom_grid).y_max = (*atom_grid).y_max + safe;
(*atom_grid).z_max = (*atom_grid).z_max + safe;
(*atom_grid).max_atom_node = (int) ((*atom_grid).mesh * (*atom_grid).mesh * 
		                  (*atom_grid).mesh / (double) min_vol_atom ) + 1;
            (*atom_grid).grid_X = 
	    (int) (((*atom_grid).x_max - (*atom_grid).x_min) / (*atom_grid).mesh ) + 1;
            
	    (*atom_grid).grid_Y = 
	    (int) (((*atom_grid).y_max - (*atom_grid).y_min) / (*atom_grid).mesh ) + 1;
            
	    (*atom_grid).grid_Z = 
	    (int) (((*atom_grid).z_max - (*atom_grid).z_min) / (*atom_grid).mesh ) + 1;
            
	    (*atom_grid).grid_size = 
	    (*atom_grid).grid_X * (*atom_grid).grid_Y * (*atom_grid).grid_Z;
//printf("%i %i %i %i\n", (*atom_grid).grid_size ,
//            (*atom_grid).grid_X , (*atom_grid).grid_Y , (*atom_grid).grid_Z);
if((*atom_grid).PBC_ON == 1)
{
(*atom_grid).dx = (*atom_grid).pbcx / (double) (*atom_grid).grid_X;
(*atom_grid).dy = (*atom_grid).pbcy / (double) (*atom_grid).grid_Y;
(*atom_grid).dz = (*atom_grid).pbcz / (double) (*atom_grid).grid_Z;
if( 
((*atom_grid).x_max - (*atom_grid).x_min) > (*atom_grid).pbcx ||
((*atom_grid).y_max - (*atom_grid).y_min) > (*atom_grid).pbcy ||
((*atom_grid).z_max - (*atom_grid).z_min) > (*atom_grid).pbcz 
)
{
printf("Warning: PBC less than system dimension...\ncontinuing....\n");
}
}
}

void atoms_on_grid(struct Atom_grid *atom_grid, struct System sistema)
{
	int i,j,k,l,n, index;
        double dx,dy,dz, safe = 0.1;
//printf("in atoms_on_grid... grid size: %i - %i %i %i\n", (*atom_grid).grid_size,
//(*atom_grid).grid_X,
//(*atom_grid).grid_Y,
//(*atom_grid).grid_Z);
for(i=0; i < (*atom_grid).grid_size; i++)
(*atom_grid).n_atom_node[i]=0;
if((*atom_grid).PBC_ON != 1)
for(n=0; n< sistema.n_atoms; n++)
{
	i = (int) ((sistema.atoms[n].coor[0] - (*atom_grid).x_min)/((*atom_grid).mesh ));
	j = (int) ((sistema.atoms[n].coor[1] - (*atom_grid).y_min)/((*atom_grid).mesh ));
	k = (int) ((sistema.atoms[n].coor[2] - (*atom_grid).z_min)/((*atom_grid).mesh ));
if(i >=  (*atom_grid).grid_X || j >=  (*atom_grid).grid_Y || k >= (*atom_grid).grid_Z || i < 0 || j < 0 || k < 0)
{
printf("%i %i %i vs %i %i %i\n", i, j, k, (*atom_grid).grid_X,(*atom_grid).grid_Y,(*atom_grid).grid_Z);
printf("%lf %lf %lf - %lf %lf %lf - %lf %lf %lf\n", 
	(*atom_grid).x_min, sistema.atoms[n].coor[0], (*atom_grid).x_max,
	(*atom_grid).y_min, sistema.atoms[n].coor[1], (*atom_grid).y_max,
	(*atom_grid).z_min, sistema.atoms[n].coor[2], (*atom_grid).z_max
	);
exit(0);
}
	index = i * (*atom_grid).grid_Y *(*atom_grid).grid_Z + j * (*atom_grid).grid_Z + k;
        (*atom_grid).atom_node[index][(*atom_grid).n_atom_node[index]] = n;
        (*atom_grid).n_atom_node[index] = (*atom_grid).n_atom_node[index] + 1;
if( (*atom_grid).n_atom_node[index] >= (*atom_grid).max_atom_node) 
{
(*atom_grid).max_atom_node =  (int) ((double) (*atom_grid).max_atom_node * sqrt(2));
for(l = 0; l< (*atom_grid).grid_size; l++)
(*atom_grid).atom_node[l] = realloc((*atom_grid).atom_node[l], (size_t) (*atom_grid).max_atom_node * sizeof(int));
}
//printf("%i - %i\n", index, (*atom_grid).n_atom_node[index]); 
//printf("%i - %i %i %i - %i\n", n, i, j, k, index);
}
else
for(n=0; n< sistema.n_atoms; n++)
{
	i = ((int) ((sistema.atoms[n].coor[0] - (*atom_grid).x_min)/(*atom_grid).dx))%(*atom_grid).grid_X ;
	j = ((int) ((sistema.atoms[n].coor[1] - (*atom_grid).y_min)/(*atom_grid).dy))%(*atom_grid).grid_Y ;
	k = ((int) ((sistema.atoms[n].coor[2] - (*atom_grid).z_min)/(*atom_grid).dz))%(*atom_grid).grid_Z ;
	index = i * (*atom_grid).grid_Y *(*atom_grid).grid_Z + j * (*atom_grid).grid_Z + k;
//printf("%i - %i %i %i - %i\n", n, i, j, k, index);
	(*atom_grid).atom_node[index][(*atom_grid).n_atom_node[index]] = n;
	(*atom_grid).n_atom_node[index]++; 
}
j=0;
//for(i=0; i< (*atom_grid).grid_size ; i++)
//if((*atom_grid).n_atom_node[i] > j) j = (*atom_grid).n_atom_node[i];
//printf("%i\n", j);

//(*atom_grid).atom_node[i] = realloc((*atom_grid).atom_node[i], (*atom_grid).n_atom_node[i] * sizeof(int)); 
for(l = 0; l< (*atom_grid).grid_size; l++)
(*atom_grid).atom_node[l] = realloc((*atom_grid).atom_node[l], (*atom_grid).n_atom_node[l] * sizeof(int));
}

/* qui non sono ancora messe le PBC */
void grid_neighbour(struct Atom_grid atom_grid, struct System sistema, 
double cutoff, struct Neighbour *neighbours, char type)
{
int i,j,k,r,s,t,index1,index2,m,n;
int i1,i2;
double x1,y1,z1,x2,y2,z2,d,tt;
double DX, DY, DZ;
for (i=0; i<sistema.n_atoms; i++) neighbours[i].n_neighbours = 0;
//printf("SONO QUI PBC = %i\n", atom_grid.PBC_ON);

for(i = 0; i< atom_grid.grid_X; i++)
for(j = 0; j< atom_grid.grid_Y; j++)
for(k = 0; k< atom_grid.grid_Z; k++)
{
index1 = i*atom_grid.grid_Y *atom_grid.grid_Z + j * atom_grid.grid_Z + k;
for(r = -1; r<= 1; r++)
for(s = -1; s<= 1; s++)
for(t = -1; t<= 1; t++)
if(atom_grid.PBC_ON != 1)
{
index2 = (i+r)*atom_grid.grid_Y *atom_grid.grid_Z + 
	 (j+s)*atom_grid.grid_Z + k+t;
	if( ((i+r) < atom_grid.grid_X ) && ((j+s) < atom_grid.grid_Y) 
		&& ((k+t) < atom_grid.grid_Z) && 
	     ((i+r) >= 0 ) && ((j+s) >= 0) && ((k+t) >= 0) )
        for(n=0; n< atom_grid.n_atom_node[index1]; n++)
	{
		i1 = atom_grid.atom_node[index1][n];
		x1 = sistema.atoms[i1].coor[0]; 
		y1 = sistema.atoms[i1].coor[1]; 
		z1 = sistema.atoms[i1].coor[2]; 
		for(m=0; m< atom_grid.n_atom_node[index2]; m++)
		{
		i2 = atom_grid.atom_node[index2][m];
		x2 = sistema.atoms[i2].coor[0]; 
		y2 = sistema.atoms[i2].coor[1]; 
		z2 = sistema.atoms[i2].coor[2];
        if(type = 'v')
	d =  (sistema.atoms[i1].radius + sistema.atoms[i2].radius + cutoff)*
             (sistema.atoms[i1].radius + sistema.atoms[i2].radius + cutoff);
        else if(type = 'd') d = cutoff*cutoff;

	if( ((x1-x2)*(x1-x2)) <= d) 
	if( ((y1-y2)*(y1-y2)) <= d) 
	if( ((z1-z2)*(z1-z2)) <= d)
        {
	tt = ((x1-x2)*(x1-x2)+(y1 - y2)*(y1 - y2)+(z1 - z2)*(z1 - z2)); 
        if(tt <= d)
        if(i1 > i2)
        {
        neighbours[i1].list[neighbours[i1].n_neighbours] = i2;
        neighbours[i2].list[neighbours[i2].n_neighbours] = i1;
        neighbours[i1].d[neighbours[i1].n_neighbours] = sqrt(tt);
        neighbours[i2].d[neighbours[i2].n_neighbours] = sqrt(tt);
        neighbours[i1].n_neighbours++;
        neighbours[i2].n_neighbours++;
        if(neighbours[i1].n_neighbours >= neighbours[i1].max_n)
         {
         neighbours[i1].max_n = neighbours[i1].max_n * 2;
         neighbours[i1].list = realloc(neighbours[i1].list, (size_t) neighbours[i1].max_n * sizeof(int));
         neighbours[i1].d = realloc(neighbours[i1].d, neighbours[i1].max_n * sizeof(double));
         }
        if(neighbours[i2].n_neighbours >= neighbours[i2].max_n)
         {
         neighbours[i2].max_n = neighbours[i2].max_n * 2;
         neighbours[i2].list = realloc(neighbours[i2].list, (size_t) neighbours[i2].max_n * sizeof(int));
         neighbours[i2].d = realloc(neighbours[i2].d, neighbours[i2].max_n * sizeof(double));
         }
	}
        }
        }
}
else if(atom_grid.PBC_ON == 1) //atom_grid.PBC_ON is set to 0 so this condition is FALSE
{
DX = atom_grid.pbcx;
DY = atom_grid.pbcy;
DZ = atom_grid.pbcz;
index2 = ((i+r+atom_grid.grid_X)%atom_grid.grid_X)*atom_grid.grid_Y *atom_grid.grid_Z + 
	 ((j+s+atom_grid.grid_Y)%atom_grid.grid_Y)*atom_grid.grid_Z + (k+t+atom_grid.grid_Z)%atom_grid.grid_Z;
if((atom_grid.grid_X > 2 || ((atom_grid.grid_X == 2) && r >= 0) || ((atom_grid.grid_X == 1) && r == 0)) 
&& (atom_grid.grid_Y > 2 || ((atom_grid.grid_Y == 2) && s >= 0) || ((atom_grid.grid_Y == 1) && s == 0))
&& (atom_grid.grid_Z > 2 || ((atom_grid.grid_Z == 2) && t >= 0) || ((atom_grid.grid_Z == 1) && t == 0))
)
        for(n=0; n< atom_grid.n_atom_node[index1]; n++)
	{
		i1 = atom_grid.atom_node[index1][n];
		x1 = sistema.atoms[i1].coor[0]; 
		y1 = sistema.atoms[i1].coor[1]; 
		z1 = sistema.atoms[i1].coor[2]; 
		for(m=0; m< atom_grid.n_atom_node[index2]; m++)
		{
		i2 = atom_grid.atom_node[index2][m];
		x2 = sistema.atoms[i2].coor[0]; 
		y2 = sistema.atoms[i2].coor[1]; 
		z2 = sistema.atoms[i2].coor[2];
        if(type = 'v')
	d =  (sistema.atoms[i1].radius + sistema.atoms[i2].radius + cutoff)*
             (sistema.atoms[i1].radius + sistema.atoms[i2].radius + cutoff);
        else if(type = 'd') d = cutoff*cutoff;
/**** SONO QUI ****/
        x2 = (fabs(x1-x2) < fabs(x1-x2+DX)) ? (x2) : (x2-DX);
        x2 = (fabs(x1-x2) < fabs(x1-x2-DX)) ? (x2) : (x2+DX);
        y2 = (fabs(y1-y2) < fabs(y1-y2+DY)) ? (y2) : (y2-DY);
        y2 = (fabs(y1-y2) < fabs(y1-y2-DY)) ? (y2) : (y2+DY);
        z2 = (fabs(z1-z2) < fabs(z1-z2+DZ)) ? (z2) : (z2-DZ);
        z2 = (fabs(z1-z2) < fabs(z1-z2-DZ)) ? (z2) : (z2+DZ);
 
	if( ((x1-x2)*(x1-x2)) <= d) 
	if( ((y1-y2)*(y1-y2)) <= d) 
	if( ((z1-z2)*(z1-z2)) <= d)
        {
	tt = ((x1-x2)*(x1-x2)+(y1 - y2)*(y1 - y2)+(z1 - z2)*(z1 - z2));
        if(tt <= d)
        if(i1 > i2)
        {
        neighbours[i1].list[neighbours[i1].n_neighbours] = i2;
        neighbours[i2].list[neighbours[i2].n_neighbours] = i1;
        neighbours[i1].d[neighbours[i1].n_neighbours] = sqrt(tt);
        neighbours[i2].d[neighbours[i2].n_neighbours] = sqrt(tt);
        neighbours[i1].n_neighbours++;
        neighbours[i2].n_neighbours++;
/*        printf("%8.3f %8.3f\n", distv(sistema.atoms[i1].coor, sistema.atoms[i2].coor ), sqrt(d)); */
        if(neighbours[i1].n_neighbours >= neighbours[i1].max_n)
         {
         neighbours[i1].max_n = neighbours[i1].max_n * 2;
         neighbours[i1].list = realloc(neighbours[i1].list, (size_t) neighbours[i1].max_n * sizeof(int));
         neighbours[i1].d = realloc(neighbours[i1].d, neighbours[i1].max_n * sizeof(double));
         }
        if(neighbours[i2].n_neighbours >= neighbours[i2].max_n)
         {
         neighbours[i2].max_n = neighbours[i2].max_n * 2;
         neighbours[i2].list = realloc(neighbours[i2].list, (size_t) neighbours[i2].max_n * sizeof(int));
         neighbours[i2].d = realloc(neighbours[i2].d, neighbours[i2].max_n * sizeof(double));
         }
        }
        }
	}
                 }
        }
}

}
}

void system_area_alloc(struct System system, struct System_area *system_area)
{
(*system_area).atoms=(struct Area *) calloc((size_t) system.n_atoms, sizeof(struct Area));
(*system_area).residues=(struct Area *) calloc((size_t) system.n_residues, sizeof(struct Area));
(*system_area).chains=(struct Area *) calloc((size_t) system.n_chains, sizeof(struct Area));
(*system_area).segments=(struct Area *) calloc((size_t) system.n_segments, sizeof(struct Area));
(*system_area).models=(struct Area *) calloc((size_t) system.n_models, sizeof(struct Area));
}

void system_area_free(struct System_area *system_area)
{
free((*system_area).atoms);
free((*system_area).residues);
free((*system_area).chains);
free((*system_area).segments);
free((*system_area).models);
}


void grid_srf_dens(
                struct System system,
                struct Atom_grid atom_grid,
                struct Flag_par flag_par,
                struct Srf *srf)
{
int i,j,k,kk,r,s,t,m,n,i1,i2,ii,nsrfpt, ir, iv, n_radii, i_srf=0;
int index1, index2, *n_srf_pt;
double acc, acc_self_chain, acc_self_res;
double phi, psi, area_sas, da = 0.0, area_ctc, d, tmpf;
double r_min, r_max, dr;
double **srfxr, **srfyr, **srfzr,*rad, dx, dist2, dist;
double *srfx, *srfy, *srfz;
double *sas_r, *sas_c;
double x1, y1, z1, x2, y2, z2, x, y, z;
double fdt,ft,fds,fs,fu;
r_max = 0.0;
r_min = DBL_MAX;
        for(i=0; i< system.n_atoms; i++)
        {
        if((system.atoms[i].radius) > r_max)
            r_max = system.atoms[i].radius ;
        if((system.atoms[i].radius ) < r_min)
            r_min = system.atoms[i].radius ;
/*        
        if((system.atoms[i].radius + flag_par.probe_radius) > r_max)
            r_max = system.atoms[i].radius + flag_par.probe_radius;
        if((system.atoms[i].radius + flag_par.probe_radius) < r_min)
            r_min = system.atoms[i].radius + flag_par.probe_radius;
*/
        }
dr = 0.01;
n_radii = (int) floor(0.000001 + (r_max - r_min)/dr) + 1;
        srfxr = (double **) calloc((size_t) n_radii, sizeof(double *));
        srfyr = (double **) calloc((size_t) n_radii, sizeof(double *));
        srfzr = (double **) calloc((size_t) n_radii, sizeof(double *));
        sas_r = (double *) calloc((size_t) n_radii, sizeof(double));
        sas_c = (double *) calloc((size_t) n_radii, sizeof(double));
        rad   = (double *) calloc((size_t) n_radii, sizeof(double));
        n_srf_pt  = (int *) calloc((size_t) n_radii, sizeof(int));
for(i = 0; i < system.n_atoms; i++) 
(*srf).at_srf[i].n_isrf  = 0;
for(i = 0; i < n_radii; i++) 
{
	n_srf_pt[i] = 4.0 * M_PI * (0.0001 + r_min + dr * (double) i ) * (0.0001 + r_min + dr * (double) i )/flag_par.sasa_area;
        srfxr[i] = (double *) calloc(n_srf_pt[i], sizeof(double));
        srfyr[i] = (double *) calloc(n_srf_pt[i], sizeof(double));
        srfzr[i] = (double *) calloc(n_srf_pt[i], sizeof(double));
}
//printf("%lf -- %lf\n", r_min, r_max); exit(0);
for(x = r_min; x<= r_max+0.0001; x = x + dr)
{
i = (int) floor((0.0001 + (x - r_min))/dr);
//printf("X %lf %i\n", x, i);
/* rad[i] = x+dr/2.0; tolta perche' spesso i raggi sono al decimo di A*/
rad[i] = x; 
fdt = M_PI * (3-sqrt(5));
ft=0.0;
fds = 2.0/(double) n_srf_pt[i];
fs = 1 - fds/2.0;
for(j = 0; j< n_srf_pt[i]; j++)
{
fu = sqrt(1-fs*fs);
srfxr[i][j] =  (x + flag_par.probe_radius ) * cos(ft)*fu;
srfyr[i][j] =  (x + flag_par.probe_radius ) * sin(ft)*fu;
srfzr[i][j] =  (x + flag_par.probe_radius ) * fs;
fs = fs -fds;
ft = ft + fdt;
}
sas_r[i] = (x + flag_par.probe_radius) * (x + flag_par.probe_radius) *4*M_PI/(double) n_srf_pt[i];
sas_c[i] = (x) * (x) *4*M_PI/(double) n_srf_pt[i];
}
for(i = 0; i< atom_grid.grid_X; i++)
for(j = 0; j< atom_grid.grid_Y; j++)
for(k = 0; k< atom_grid.grid_Z; k++)
{
     index1 = i*atom_grid.grid_Y *atom_grid.grid_Z + j * atom_grid.grid_Z + k;
        for(n=0; n< atom_grid.n_atom_node[index1]; n++)
          {
                i1 = atom_grid.atom_node[index1][n];
                if ((system.atoms[i1].radius >= r_min) && (system.atoms[i1].radius <= r_max)) 
                {
                ir = (int) floor( (system.atoms[i1].radius - r_min)/dr + 0.0001);
                area_sas = sas_r[ir];
//                printf("%i %i\n",ir, n_radii); 
//exit(0);
                srfx = srfxr[ir]; srfy = srfyr[ir]; srfz = srfzr[ir];
                }
                else printf("I don't know ATOM: element %s: %s in res %s \n", system.atoms[i1].element,  system.atoms[i1].at_name, system.atoms[i1].res_name);
/*        if(strcmp(system.atoms[i1].element,"UN")) */
        if(system.atoms[i1].radius < r_min || system.atoms[i1].radius > r_max)
            {
            printf("%s %s %i %c has radius %11.3f... reset to 1.5\n",
            system.atoms[i1].at_name,
            system.atoms[i1].res_name,
            system.atoms[i1].res_n,
            system.atoms[i1].chain,
            system.atoms[i1].radius);
            system.atoms[i1].radius = 1.5;
            }
        for(ii = 0; ii<n_srf_pt[ir]; ii++)
        {
        x1 = system.atoms[i1].coor[0] + srfx[ii];
        y1 = system.atoms[i1].coor[1] + srfy[ii];
        z1 = system.atoms[i1].coor[2] + srfz[ii];
//	printf("srfpt: %lf %lf %lf of atom %i ", x1, y1, z1, i1);
        acc_self_res = acc_self_chain = acc = 1.0;
        for(r = -1; r<= 1; r++)
        for(s = -1; s<= 1; s++)
        for(t = -1; t<= 1; t++)
        {
        index2 = (i+r)*atom_grid.grid_Y * atom_grid.grid_Z + (j+s)*atom_grid.grid_Z + k+t;
        if(
           ((i+r) < atom_grid.grid_X ) && ((j+s) < atom_grid.grid_Y) &&
           ((k+t) < atom_grid.grid_Z) && ((i+r) >= 0 ) && ((j+s) >= 0) &&
           ((k+t) >= 0)
          )
        for(m=0; m< atom_grid.n_atom_node[index2]; m++)
        if(atom_grid.atom_node[index2][m] != i1)
        {
        i2 = atom_grid.atom_node[index2][m];
//        printf("%i %i\n", i1, i2); 
        x2 = system.atoms[i2].coor[0];
        y2 = system.atoms[i2].coor[1];
        z2 = system.atoms[i2].coor[2];
        d = system.atoms[i2].radius+flag_par.probe_radius;
        x = x1-x2;
        if(fabs(x) < d)
         {
          y = y1-y2;
           if(fabs(y) < d)
            {
             z = z1-z2;
              if(fabs(z) < d)
               if(i1!=i2)
               {
               dist2 = x*x + y*y + z*z;
               if(dist2 < (d*d))
                  {
                  acc = 0;
//	printf("excluded by: %lf %lf %lf of atom %i\n", x2, y2, z2, i2);
                  if(
                     ((system.atoms[i2].chain == system.atoms[i1].chain) && 
                      !(strcmp(system.atoms[i1].segid,system.atoms[i2].segid)))
                    )
                    {
                    acc_self_chain = 0; 
                    if(system.atoms[i2].res_n == system.atoms[i1].res_n)
                    if(system.atoms[i2].res_ins == system.atoms[i1].res_ins)
                    acc_self_res = 0; 
                    }
                  }
              }
            }
          }
        if(acc == 0) // break the for cycle
        m=atom_grid.n_atom_node[index2];
         }
	}
        if(acc > 0.0)
        {
        if(i_srf >= ((*srf).alloc_srf_pt)) {
printf("Attempting to write more surface points than allocated, I will reallocate surface points...\n");
         (*srf).alloc_srf_pt = (int) ((double) (*srf).alloc_srf_pt * sqrt(2));
         (*srf).srf_pt = realloc((*srf).srf_pt, sizeof(struct Srf_pt) * (*srf).alloc_srf_pt);
         for(kk=i_srf; kk< (*srf).alloc_srf_pt ; kk++)
            {
            (*srf).srf_pt[kk].r = (double *) calloc((size_t) 3, sizeof(double));
            (*srf).srf_pt[kk].vec = (double *) calloc((size_t) 3, sizeof(double));
            }
         if((*srf).srf_pt == NULL) printf("Could not allocate %i srf points\n...exiting\n", (*srf).alloc_srf_pt ); exit(0);
         }
        (*srf).srf_pt[i_srf].a = area_sas  * acc ;
        (*srf).srf_pt[i_srf].r[0] = x1;
        (*srf).srf_pt[i_srf].r[1] = y1;
        (*srf).srf_pt[i_srf].r[2] = z1;
        (*srf).srf_pt[i_srf].vec[0] = srfx[ii]/(rad[ir] + flag_par.probe_radius);
        (*srf).srf_pt[i_srf].vec[1] = srfy[ii]/(rad[ir] + flag_par.probe_radius);
        (*srf).srf_pt[i_srf].vec[2] = srfz[ii]/(rad[ir] + flag_par.probe_radius);
        (*srf).srf_pt[i_srf].at = i1;
        (*srf).at_srf[i1].isrf[ (*srf).at_srf[i1].n_isrf++] = i_srf;
        i_srf++;
        }
   }
  }
}
for(kk = i_srf; kk<(*srf).alloc_srf_pt; kk++)
{
free((*srf).srf_pt[kk].r);
free((*srf).srf_pt[kk].vec);
}
(*srf).n_srf_pt_tot = i_srf;
(*srf).alloc_srf_pt = i_srf;
         (*srf).srf_pt = (struct Srf_pt *) realloc((*srf).srf_pt, (size_t) (sizeof(struct Srf_pt) * (*srf).alloc_srf_pt));
//for(i = 0; i < n_radii; i++) 
//printf("%lf %i\n", r_min + dr * (double) i,  n_srf_pt[i]);
for(i = 0; i < n_radii; i++) 
{
        free(srfxr[i]); 
        free(srfyr[i]); 
        free(srfzr[i]); 
}
        free(srfxr); 
        free(srfyr); 
        free(srfzr); 
        free(sas_r); 
//for(i = 0; i < system.n_atoms; i++) 
//printf("srf_at: %i -- %i\n", i, (*srf).at_srf[i].n_isrf);
}

void sort_neighbours(struct Neighbour *neighbours, int beg, int end)
{
   int  l,r,p;
   double pivd, tmpd;
   int pivi, tmpi;

    while (beg<end)    // This while loop will substitude the second recursive call
    {
        l = beg; p = (beg+end)/2; r = end;

        pivd = (*neighbours).d[p];
        pivi = (*neighbours).list[p];

        while (1)
        {
          // while ((l<=r) && ((*neighbours).list[l] <= pivi || ((*neighbours).list[l] == pivi && (*neighbours).d[l] <= pivd ))) l++;
         // while ((l<=r) && ((*neighbours).list[r] > pivi || ((*neighbours).list[r] == pivi && (*neighbours).d[r] > pivd ))) r--;
            while ((l<=r) && ((*neighbours).d[l] <= pivd )) l++;
            while ((l<=r) && ((*neighbours).d[r] > pivd )) r--;

            if (l>r) break;

            tmpd=(*neighbours).d[l]; (*neighbours).d[l]=(*neighbours).d[r]; (*neighbours).d[r]=tmpd;
            tmpi=(*neighbours).list[l]; (*neighbours).list[l]=(*neighbours).list[r]; (*neighbours).list[r]=tmpi;

            if (p==r) p=l;
            
            l++; r--;
        }

        (*neighbours).d[p]=(*neighbours).d[r]; (*neighbours).d[r]=pivd;
        (*neighbours).list[p]=(*neighbours).list[r]; (*neighbours).list[r]=pivi;
        r--;

        // Select the shorter side & call recursion. Modify input param. for loop
        if ((r-beg)<(end-l))   
        {
            sort_neighbours(neighbours, beg, r);
            beg=l;
        }
        else
        {
            sort_neighbours(neighbours, l, end);
            end=r;
        }
    }   
}

void sort_neighbours_by_dist(struct Neighbour *neighbours, int beg, int end)
{
   int  l,r,p;
   double pivd, tmpd;
   int pivi, tmpi;

    while (beg<end)    // This while loop will substitude the second recursive call
    {
        l = beg; p = (beg+end)/2; r = end;

        pivd = (*neighbours).d[p];
        pivi = (*neighbours).list[p];

        while (1)
        {
            while ((l<=r) && ((*neighbours).d[l] <= pivd )) l++;
            while ((l<=r) && ((*neighbours).d[r] > pivd )) r--;

            if (l>r) break;

            tmpd=(*neighbours).d[l]; (*neighbours).d[l]=(*neighbours).d[r]; (*neighbours).d[r]=tmpd;
            tmpi=(*neighbours).list[l]; (*neighbours).list[l]=(*neighbours).list[r]; (*neighbours).list[r]=tmpi;

            if (p==r) p=l;
            
            l++; r--;
        }

        (*neighbours).d[p]=(*neighbours).d[r]; (*neighbours).d[r]=pivd;
        (*neighbours).list[p]=(*neighbours).list[r]; (*neighbours).list[r]=pivi;
        r--;

        // Select the shorter side & call recursion. Modify input param. for loop
        if ((r-beg)<(end-l))   
        {
            sort_neighbours_by_dist(neighbours, beg, r);
            beg=l;
        }
        else
        {
            sort_neighbours_by_dist(neighbours, l, end);
            end=r;
        }
    }   
}

void print_grid_info(struct Atom_grid atom_grid)
{
printf("grid (X,Y,Z): %i %i %i -- grid size: %i -- grid mesh: %lf\n", atom_grid.grid_X, atom_grid.grid_Y, atom_grid.grid_Z, atom_grid.grid_size, atom_grid.mesh);
if(atom_grid.PBC_ON)
printf("x: %lf -- %lf    dx: %lf\ny: %lf -- %lf    dy: %lf\nz: %lf -- %lf    dz: %lf\n",
atom_grid.x_min, atom_grid.x_max, atom_grid.dx,
atom_grid.y_min, atom_grid.y_max, atom_grid.dy,
atom_grid.z_min, atom_grid.z_max, atom_grid.dz);
else 
printf("x: %lf -- %lf\ny: %lf -- %lf\nz: %lf -- %lf\n",
atom_grid.x_min, atom_grid.x_max,
atom_grid.y_min, atom_grid.y_max,
atom_grid.z_min, atom_grid.z_max);
}

void srf_free(struct Srf *srf, int n_atoms)
{
int i;
for(i=0;i<(*srf).alloc_srf_pt;i++)
{
free((*srf).srf_pt[i].r);
free((*srf).srf_pt[i].vec);
} 
free((*srf).srf_pt);
for(i=0;i<n_atoms;i++)
free((*srf).at_srf[i].isrf);
free((*srf).at_srf);
}
