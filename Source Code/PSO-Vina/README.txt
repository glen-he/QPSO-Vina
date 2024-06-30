=======================================================================
PSOVina 2.0
Giotto H. K. Tai & Shirley W. I. SIU                                   
                                                                       
Computational Biology and Bioinformatics Lab (CBBio)                        
University of Macau                                                    
                                                                       
http://cbbio.cis.umac.mo/software/psovina
=======================================================================

0. WHAT IS IT?
   A fast protein-ligand docking tool based on the implementation of 
   AutoDock Vina and Particle Swarm Optimization (PSO) method. This version
   came along from our initial psovina (psovina-1.1) in 2015, then added 
   2-stage local search to shorten the expensive local search (psovina2ls-1.2)
   in 2016. 

   In this version, we further enhance pose prediction by embedding the
   Singer chaotic function into the global PSO search to achieve a 5 to 6-fold
   speedup. This makes PSOVina 2.0 an ideal tool for virtual screening task
   or other parameter optimization purpose, whenever many iterative docking
   experiments are needed to be done. 

   - PSOVina 2.0: Chaos-embedded with enhanced pose prediction with speed
   - PSOVina 2LS: Designed 2-stage local search, the fastest docking tool
   - PSOVina 1.0: Hybrid AutoDock Vina and PSO 

1. INSTALLATION
   For successful compilation, please install Boost (version 1.59.0) 
   from http://www.boost.org. For preparing molecules for docking, 
   please install AutoDockTools (ADT) from http://mgltools.scripps.edu.

   The installation basically follows the installation of AutoDock Vina. 
   The steps are simple:

     a. unpack the files
     b. cd psovina-2.0/build/<your-platform>/release
     c. modify Makefile to suit your system setting
     d. type "make" to compile
    
    The binary psovina will be generated at the current directory. You can 
    copy this binary to a directory in your PATH e.g. /usr/local/bin, or add
    the path of the current directory to your PATH.

2. RUNNING PSOVINA

   You can run psovina as the way you run vina but additional three 
   parameters (optional) are used to specify how the PSO algorithm performs 
   searching:

   % <path-to-psovina>/psovina

   PSO parameters (optional):
       --num_particles arg (=8)      Number of particles
       --w arg (=0.36)               Inertia weight
       --c1 arg (=0.99)              Cognitive weight 
       --c2 arg (=0.99)              Social weight 

   For example, docking Kifunensine in the Mannosidase enzyme (PDBID 1ps3 from
   the PDBbind v2012 dataset) using PSOVina with default PSO parameters in a 
   8-core computer and return the lowest energy prediction:

   % <path-to-AutoDockTools>/prepare_ligand4.py -l 1ps3_ligand.mol2 \
     -o 1ps3_ligand.pdbqt -A 'hydrogens' -U 'nphs_lps_waters'

   % <path-to-AutoDockTools>/prepare_receptor4.py -r 1ps3_protein.pdb \
     -o 1ps3_protein.pdbqt -A 'hydrogens' -U 'nphs_lps_waters' 

   % <path-to-psovina>/psovina \
     --receptor 1ps3_protein.pdbqt --ligand 1ps3_ligand.pdbqt \
     --center_x  31.951 --center_y 65.5053 --center_z 7.63888 \
     --size_x    33.452 --size_y   27.612  --size_z   35.136  \ 
     --cpu 8  
   

3. DEVELOP PSOVINA

   If you are interested in the source code of PSOVina for any academic
   purposes, please note that the following files were newly developed   
   in our work or modified based on PSOVina:
	src/lib/pso.h
	src/lib/pso.cpp
	src/lib/pso_mutate.cpp
	src/lib/pso_search.cpp
	src/lib/quasi_newton.h
	src/lib/quasi_newton.cpp

4. CITATION
   Please cite our paper if you have used any version of PSOVina for your work. 

   Hio Kuan Tai, Siti Azma Jusoh, and Shirley W. I. Siu  
   PSOVina2.0: Improving Protein-Ligand Docking by Chaos-Embedded Particle Swarm Optimization
   (Submitted) 

   Hio Kuan Tai, Hin Lin and Shirley W. I. Siu 
   Improving the efficiency of PSOVina for protein-ligand docking by two-stage local search
   2016 IEEE Congress on Evolutionary Computation (CEC), Vancouver, BC, 2016, pp. 770-777.

   Marcus C. K. Ng, Simon Fong, and Shirley W. I. Siu
   PSOVina: The Hybrid Particle Swarm Optimization Algorithm for Protein-Ligand Docking
   Journal of Bioinformatics and Computational Biology (JBCB) 13 (3), 1541007, 2015.
   
   Please check out our homepage for the updated citation.

5. CONTACT US
    
   Developers: 
     PSOVina 2.0: Giotto H. K. Tai (giottotai@yahoo.com.hk)
     PSOVina 2LS: Allan H. Lin     (lhang33@126.com)
     PSOVina 1.0: Marcus C. K. Ng  (marcus.ckng@gmail.com)

   Project P.I.: 
     Shirley W. I. Siu (shirleysiu@umac.mo)

   Computational Biology and Bioinformatics Lab, University of Macau
   http://cbbio.cis.umac.mo
   http://www.cis.umac.mo/~shirleysiu
