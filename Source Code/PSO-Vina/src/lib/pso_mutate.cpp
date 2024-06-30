/*
       PSOVina version 2.0  

        This file is revised from monte_carlo.cpp in AutoDock Vina.

        Authors: Giotto H. K. TAI  <giottotai@yahoo.com.hk>

                 Shirley W. I. Siu <shirleysiu@umac.mo>


        Computational Biology and Bioinformatics Lab

        University of Macau

        http://cbbio.cis.umac.mo

*/
/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

#include "pso_mutate.h"
#include <math.h>
#include <time.h>
#define PI 3.14159265

sz pso_count_mutable_entities(const conf& c) {
	sz counter = 0;
	VINA_FOR_IN(i, c.ligands)
		counter += 2 + c.ligands[i].torsions.size();
	VINA_FOR_IN(i, c.flex)
		counter += c.flex[i].torsions.size();
	return counter;
}

void pso_mutate_conf(output_type& candidate, output_type& candidate_1, const model& m, fl amplitude, rng& generator, pso* particle, double* PersonalBest, const precalculate& p ,const igrid& ig,change& g,const vec& v,quasi_newton& quasi_newton_par,int step) { // ONE OF: 2A for position, similar amp for orientation, randomize torsion
	int shrink = 1; //rough factor:1/10 = 0.1
	output_type tmp_2 = candidate;
	sz mutable_entities_num = pso_count_mutable_entities(candidate.c);
	if(mutable_entities_num == 0) return;
	int which_int = random_int(0, int(mutable_entities_num - 1), generator);
	VINA_CHECK(which_int >= 0);
	sz which = sz(which_int);
	VINA_CHECK(which < mutable_entities_num);
	
	srand((unsigned)time(NULL));
	double r_cm;
	r_cm = (double) rand() / (RAND_MAX + 1.0);
	
	double w_cm;
	w_cm = (double) rand() / (RAND_MAX + 1.0);

	int y;
	VINA_FOR_IN(i, candidate.c.ligands) {
			
		model tmp_m = m;
		const vec authentic_v(1000, 1000, 1000);
		//Take part the position (either position or orientation or torsion)
		if(which == 0) {
			for (y=0;y<particle->number;y++)
			{
				candidate.c.ligands[i].rigid.position = particle->getCurrentPosition(y);
				candidate.c.ligands[i].rigid.orientation = particle->getCurrentOrientation(y);
				for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
					candidate.c.ligands[i].torsions[z] = particle->getCurrentTorsion(y,z);
				
				//rough local search
				quasi_newton_par(tmp_m, p, ig, candidate, g, v, shrink);
				
				//set the personal best(energy value and position);
				if(candidate.e < particle->getPersonalBest(y) || step <= 18) //Cr = 18
				{
					if(candidate.e < particle->getPersonalBest(y))
						particle->updatePersonalBest(y,candidate.e);
					tmp_2 = candidate;
					
					//full local search
					quasi_newton_par(tmp_m, p, ig, tmp_2, g, v);
					
					if (tmp_2.e < PersonalBest[y])
					{
						PersonalBest[y] = tmp_2.e;
						particle->updateBestPosition(y,tmp_2.c.ligands[i].rigid.position);
						particle->updateBestOrientation(y,tmp_2.c.ligands[i].rigid.orientation);
						for(int z=0;z<tmp_2.c.ligands[i].torsions.size();z++)
							particle->updateBestTorsion(y, tmp_2.c.ligands[i].torsions[z],z);
					}
						
					//set the global best(energy value and position);
					if(tmp_2.e < pso::gbest_fit)
					{
						particle->updateGlobalBest_1(tmp_2.e);
						
						pso::gbest_position = tmp_2.c.ligands[i].rigid.position;
						pso::gbest_orientation = tmp_2.c.ligands[i].rigid.orientation;
						for(int z=0;z<tmp_2.c.ligands[i].torsions.size();z++)
							pso::gbest_torsion[z] = tmp_2.c.ligands[i].torsions[z];
					}
				}
				
				//update chaotic map
				r_cm = 1.07*(7.86*r_cm-23.31*r_cm*r_cm+28.75*r_cm*r_cm*r_cm-13.302875*r_cm*r_cm*r_cm*r_cm);
				w_cm = 1.07*(7.86*w_cm-23.31*w_cm*w_cm+28.75*w_cm*w_cm*w_cm-13.302875*w_cm*w_cm*w_cm*w_cm);
				
				//update each particle in every dimension
				particle->updateVelocity(generator,y,r_cm,w_cm);
				
				//compute the new position;
				particle->computeNewPositions(y);
			}

			for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
				candidate_1.c.ligands[i].torsions[z] = pso::gbest_torsion[z];
			candidate_1.c.ligands[i].rigid.orientation = pso::gbest_orientation;
			candidate_1.c.ligands[i].rigid.position = pso::gbest_position;
			candidate_1.e = pso::gbest_fit;
			return; 
			}
			
			
			--which;
			//Take part orientation (either position or orientation or torsion)
			if(which == 0) { 
				fl gr = m.gyration_radius(i); 
				if(gr > epsilon_fl) { // FIXME? just doing nothing for 0-radius molecules. do some other mutation?
			
					for (y=0;y<particle->number;y++)
					{
						candidate.c.ligands[i].rigid.position = particle->getCurrentPosition(y);
						candidate.c.ligands[i].rigid.orientation = particle->getCurrentOrientation(y);
						for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
							candidate.c.ligands[i].torsions[z] = particle->getCurrentTorsion(y,z);

						//rough local search
						quasi_newton_par(tmp_m, p, ig, candidate, g, v, shrink);

					
						//set the personal best(energy value and position);
						if(candidate.e < particle->getPersonalBest(y) || step <= 18)
						{
							if(candidate.e < particle->getPersonalBest(y))
								particle->updatePersonalBest(y,candidate.e);
							tmp_2 = candidate;

							//full local search
							quasi_newton_par(tmp_m, p, ig, tmp_2, g, v);
							
							if(tmp_2.e < PersonalBest[y])
							{
								PersonalBest[y] = tmp_2.e;
								particle->updateBestPosition(y,tmp_2.c.ligands[i].rigid.position);
								particle->updateBestOrientation(y,tmp_2.c.ligands[i].rigid.orientation);
								for(int z=0;z<tmp_2.c.ligands[i].torsions.size();z++)
									particle->updateBestTorsion(y, tmp_2.c.ligands[i].torsions[z],z);
							}

							//set the global best(energy value and position);
							if(tmp_2.e< pso::gbest_fit)
							{
								particle->updateGlobalBest_1(tmp_2.e);						
								pso::gbest_position = tmp_2.c.ligands[i].rigid.position;
								pso::gbest_orientation = tmp_2.c.ligands[i].rigid.orientation;
								for(int z=0;z<tmp_2.c.ligands[i].torsions.size();z++)
									pso::gbest_torsion[z] = tmp_2.c.ligands[i].torsions[z];
							}
						}

					r_cm = 1.07*(7.86*r_cm-23.31*r_cm*r_cm+28.75*r_cm*r_cm*r_cm-13.302875*r_cm*r_cm*r_cm*r_cm);
					w_cm = 1.07*(7.86*w_cm-23.31*w_cm*w_cm+28.75*w_cm*w_cm*w_cm-13.302875*w_cm*w_cm*w_cm*w_cm);
					//update each particle in every dimension
					particle->updateVelocityO(generator,y,r_cm,w_cm);
					//compute the new position;
					particle->computeNewOrientation(y);
					}
				}

				for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
					candidate_1.c.ligands[i].torsions[z] = pso::gbest_torsion[z];
				candidate_1.c.ligands[i].rigid.orientation = pso::gbest_orientation;
				candidate_1.c.ligands[i].rigid.position = pso::gbest_position;
				candidate_1.e = pso::gbest_fit;
				return; 
			}
			
			/*Torsions*/
			--which;
			if(which < candidate.c.ligands[i].torsions.size()) {
				for (y=0;y<particle->number;y++)
				{
					candidate.c.ligands[i].rigid.position = particle->getCurrentPosition(y);
					candidate.c.ligands[i].rigid.orientation = particle->getCurrentOrientation(y);
					for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
						candidate.c.ligands[i].torsions[z] = particle->getCurrentTorsion(y,z);
					
					//rough local search
					quasi_newton_par(tmp_m, p, ig, candidate, g, v, shrink);

					//set the personal best(energy value and position);
					if(candidate.e < particle->getPersonalBest(y) || step <= 18)
					{
						if(candidate.e < particle->getPersonalBest(y))
							particle->updatePersonalBest(y,candidate.e);
						tmp_2 = candidate;
	
						//full local search
						quasi_newton_par(tmp_m, p, ig, tmp_2, g, v);
						
						if(tmp_2.e < PersonalBest[y])
						{
							PersonalBest[y] = tmp_2.e;
							particle->updateBestPosition(y,tmp_2.c.ligands[i].rigid.position);
							particle->updateBestOrientation(y,tmp_2.c.ligands[i].rigid.orientation);
							for(int z=0;z<tmp_2.c.ligands[i].torsions.size();z++)
								particle->updateBestTorsion(y, tmp_2.c.ligands[i].torsions[z],z);
						}
						
						
						//set the global best(energy value and position);
						if(tmp_2.e < pso::gbest_fit)
						{
							particle->updateGlobalBest_1(tmp_2.e);
							pso::gbest_position = tmp_2.c.ligands[i].rigid.position;
							pso::gbest_orientation = tmp_2.c.ligands[i].rigid.orientation;
							for(int z=0;z<tmp_2.c.ligands[i].torsions.size();z++)
								pso::gbest_torsion[z] = tmp_2.c.ligands[i].torsions[z];
						}
					}
					
				//update chaotic map
				r_cm = 1.07*(7.86*r_cm-23.31*r_cm*r_cm+28.75*r_cm*r_cm*r_cm-13.302875*r_cm*r_cm*r_cm*r_cm);
				w_cm = 1.07*(7.86*w_cm-23.31*w_cm*w_cm+28.75*w_cm*w_cm*w_cm-13.302875*w_cm*w_cm*w_cm*w_cm);
				
				//update each particle in every dimension
				particle->updateVelocityT(generator,y,which,r_cm,w_cm);
				
				//compute the new position;
				particle->computeNewTorsion(y,generator,which);
				}
				
				for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
					candidate_1.c.ligands[i].torsions[z] = pso::gbest_torsion[z];
				candidate_1.c.ligands[i].rigid.orientation = pso::gbest_orientation;
				candidate_1.c.ligands[i].rigid.position = pso::gbest_position;
				candidate_1.e = pso::gbest_fit;
				return; 
			}
			which -= candidate.c.ligands[i].torsions.size();
		}
		
	
	VINA_FOR_IN(i, candidate.c.flex) {
		if(which < candidate.c.flex[i].torsions.size()) { candidate.c.flex[i].torsions[which] = random_fl(-pi, pi, generator); return; }
		which -= candidate.c.flex[i].torsions.size();
	}
}
