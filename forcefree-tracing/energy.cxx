    // safe allocate
#define ALLOCATE(A,LEN,TYPE)                                             \
  if ( !((A)=(TYPE *)malloc((size_t)(LEN)*sizeof(TYPE))) ) ERROR(("Cannot allocate."));

  {  // start energy diagnostics
    
    static int energy_first_time_called=1;
    
    FileIO fedump,f;                                      // file to dump energy data
    FileIOStatus status;                                // signal from instantiations of FileIO class
    
    static float *dist;                                 // array to hold the distribution function

    static float *edist;                                 // array to hold the distribution function
    
    static int nex;                                 // number of energy bands
    static double emax;                                 // max energy (in units of vthe^2) 


    static std::vector<edata *> edParams;               // vector to hold species parameters, such as vth

    static int64_t total_cells;

    int64_t array_length, icell, i;

    static int nsp;
    particle_t *p;
    species_t *sp; 
    int k, isp; 

    char fname [256];

    float * wp;

	  double vth, dke;

		int nbin = 500;

		float eminp = 0.001;
		float emaxp = 100;
    float dloge;
		float de;

  if (should_dump(ehydro)){


	if (energy_first_time_called ) {

		sim_log("initializing the energy diagnostics"); 
		
		nex = global->nex;
		emax = global->emax;

		total_cells=(grid->nx+2)*(grid->ny+2)*(grid->nz+2);    

		array_length = nex*total_cells;

		ALLOCATE(dist, array_length, float);
		
		ALLOCATE(edist, nbin, float);

		dloge = (log10(emaxp)-log10(eminp))/(nbin);

		edParams.push_back(&global->ede);
		edParams.push_back(&global->edi);

		nsp = 2;
					
		energy_first_time_called = 0;
		
	} // first time called

	
	for (isp=0; isp<nsp; isp++) {   // loop over species
		
		sp = find_species_id( edParams.at(isp)->sp_id, species_list);
		sim_log("computing the distribution function for species "<< sp->name); 
		vth = edParams.at(isp)->vth;
		dke = emax*(vth*vth/2.0)/nex;
		
		p = sp->p;   // header of the particle array
		
		for (icell=0; icell< total_cells; icell++)
		for ( k=0; k<nex; k++) dist[k*total_cells+icell]=0; //zero out distribution function
		
		for ( k=0; k<nbin; k++) edist[k]=0; //zero out distribution function

		for ( i=0; i<sp->np; ++i ) {                 // loop over particles

			double gam2 = 1.0+p->ux*p->ux + p->uy*p->uy + p->uz*p->uz;
			double ke = sqrt(gam2) - 1.0;

			k = int(ke/dke);                          // energy bin
			//if (k > (nex-1)) k = nex-1;               // everything with energy > emax goes into the last bin
			if (k<=nex-1 && k>=0) dist[k*total_cells+p->i]++;               // incerement the corresponding bin for cell p->i

			de = ke*dloge;
			k = (log10(ke)-log10(eminp))/dloge + 1;
			if (k<nbin-1 && k>=0) edist[k]++;

			p++;                                      // next particle

		} 


		// normalize the distribution function 

		//{ // start normalization
		
		//int ix,iy,iz, ixn, iyn, izn , gcell;
		//int nxp2 = grid->nx + 2;
		//int nyp2 = grid->ny + 2;
		//int nzp2 = grid->nz + 2;

		//for (iz=0; iz<nzp2; iz++)
    //for (iy=0; iy<nyp2; iy++)
	  //for (ix=0; ix<nxp2; ix++) {
		
		//np=0;  // particle counter 
		
		//icell = LOCAL_CELL_ID(ix,iy,iz);
		
		//for ( k=0; k<nex; k++) np += dist[k*total_cells+icell];  // count the particles in the cell
		//if (np>0) for ( k=0; k<nex; k++) dist[k*total_cells+icell] /= np;  // normalize the distribution 
		
		
		// is this a ghost cell ?

		//gcell = (ix==0) || (ix==nxp2-1) || (iy==0) || (iy==nyp2-1) || (iz==0) || (iz==nzp2-1) ;
		
		// if yes, assign the value from a neighboring cell
		
		//if (gcell) {
			
			// find the neighboring cell
 
	//		ixn=ix; iyn=iy; izn=iz;
	 
		//	if (ix==0) ixn++;
	//		if (ix==nxp2-1) ixn--;
		//	if (iy==0) iyn++;
	//		if (iy==nyp2-1) iyn--;
		//	if (iz==0) izn++;
	//		if (iz==nzp2-1) izn--;
			
		//	int64_t nid = LOCAL_CELL_ID (ixn,iyn,izn);
			
		//	for ( k=0; k<nex; k++) dist[k*total_cells+icell] = dist[k*total_cells+nid];
			
		//}

	//}


		//} // end normalization
		
		// dump the file
		
		sim_log(" writing the distribution function to file "); 	 
		sprintf(fname,HYDRO_FILE_FORMAT,step,edParams.at(isp)->fname,step,(int)rank()); 
		sim_log("appendig data to "<<fname);

		status = fedump.open( fname, io_append );         // open the file 
		if ( status==fail ) ERROR(("Could not open file."));   // error check
		
		// dump all the energy bands

		for ( k=0; k<nex; k++) {
      wp = dist+(k*total_cells);
      fedump.write( wp, total_cells );  
		}	

		fedump.close();

    sprintf(fname,SPEC_FILE_FORMAT,step,edParams.at(isp)->fname,step,(int)rank());
      status = f.open(fname,io_write);
      if ( status==fail ) ERROR(("Could not open spectrum file."));   // error check
      f.write( edist, nbin );
      f.close();
     
	
	
    } // end species loop
    
    } // if (should_dump(hydro))

  }     // end energy diagnostics

