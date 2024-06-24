/*
        PSOVina version 2.0 

        Authors: Giotto H. K. TAI  <giottotai@yahoo.com.hk>

                 Shirley W. I. Siu <shirleysiu@umac.mo>


        Computational Biology and Bioinformatics Lab

        University of Macau

        http://cbbio.cis.umac.mo

*/
#include "pso.h"
#include "random.h"

    // the vector of degree of freedom
	qt pso::gbest_orientation;
	fl* pso::gbest_torsion;
	vec pso::gbest_position;

	double pso::gbest_fit;

	pso::pso(int num_birds,double w,double c1,double c2,const vec corner1,const vec corner2, rng& g,conf& c)
	{
		  sz torsionSize = c.ligands[0].torsions.size();
		  this->w = w;
		  this->c1 = c1;
		  this->c2 = c2;
		  this->number = num_birds;
		  this->g = g;
		  this->corner1[0] = corner1[0];			//minmum
		  this->corner1[1] = corner1[1];
		  this->corner1[2] = corner1[2];
		  this->corner2[0] = corner2[0];			//maximum
		  this->corner2[1] = corner2[1];
		  this->corner2[2] = corner2[2];
		  this->torsionSize = (int)torsionSize;
		  this->R1Max_ = 1;
		  this->R1Min_ = 0;
		  this->R2Max_ = 1;
		  this->R2Min_ = 0;
		  pso::gbest_torsion = new fl[torsionSize];
		  init(g,c);
	}
	
	void pso::init(rng &g,conf& c)
	{

		int i;
		for(i=0;i<this->number;i++)
		{
			bird single_bird;

			single_bird.pbest_fit = 1.7976931348623158e+308;
			single_bird.tmp_fit = 1.7976931348623158e+308;
			
			//set position part
			single_bird.velocity = random_in_box(this->corner1,this->corner2,g);
			single_bird.current_position = random_in_box(this->corner1,this->corner2,g);
			
			//set orientation part
			single_bird.vO = random_inside_sphere(g);
			qt tmp_o = c.ligands[0].rigid.orientation;
			quaternion_increment(tmp_o,  random_inside_sphere(g));
			single_bird.current_orientation = tmp_o;
			
			//init. the array for the number of torsion
			single_bird.current_torsion=new fl[this->torsionSize];
			single_bird.vT=new fl[this->torsionSize];
			single_bird.pbest_torsion=new fl[this->torsionSize];
			
			for(int x=0;x<this->torsionSize;x++)						//init. all the torsion that the ligand has
			{
				single_bird.vT[x] = random_fl(-pi, pi, g);
				single_bird.current_torsion[x] = random_fl(-pi, pi, g);
			
			}
			
			particle.push_back(single_bird);					
		}
		
		
		pso::gbest_fit = 1.7976931348623158e+308;

	}
	
	
	void pso::updateVelocity(rng& generator,int i,double cm,double l)
	{
			particle[i].velocity[0] = particle[i].velocity[0]*l+c1*cm*(particle[i].pbest_pos[0]-particle[i].current_position[0])+c2*(1-cm)*(pso::gbest_position[0]-particle[i].current_position[0]);
			particle[i].velocity[1] = particle[i].velocity[1]*l+c1*cm*(particle[i].pbest_pos[1]-particle[i].current_position[1])+c2*(1-cm)*(pso::gbest_position[1]-particle[i].current_position[1]);
			particle[i].velocity[2] = particle[i].velocity[2]*l+c1*cm*(particle[i].pbest_pos[2]-particle[i].current_position[2])+c2*(1-cm)*(pso::gbest_position[2]-particle[i].current_position[2]);
				
	}
	
	void pso::updateVelocityO(rng& generator,int i,double cm,double l)
	{
	
		    qt p1 = particle[i].pbest_orientation-particle[i].current_orientation;
		    qt p2 = pso::gbest_orientation-particle[i].current_orientation;
			particle[i].vO = particle[i].vO*l+c1*cm*quaternion_to_angle(p1)+c2*(1-cm)*quaternion_to_angle(p2);
			
	}
	
	void pso::updateVelocityT(rng& generator,int i,sz which,double cm,double l)
	{
			particle[i].vT[which] = particle[i].vT[which]*l+c1*cm*(particle[i].pbest_torsion[which]-particle[i].current_torsion[which])+c2*(1-cm)*(pso::gbest_torsion[which]-particle[i].current_torsion[which]);
	}
	
	void pso::computeNewPositions(int i, int max_iter, int current_iter, rng& generator)
	{
		double beta = (1 - 0.5) * (max_iter - current_iter) / static_cast<fl>(max_iter) + 0.5;

		int num = pso::number;
		vec mbest_tmp;
		mbest_tmp.assign(0.0);
		for (int j = 1; j <= num; ++j)
		{
			mbest_tmp += particle[j].pbest_pos;
		}
		vec mbest = mbest_tmp / static_cast<double>(num);
	
		double phi = random_double(0, 1, g);
		vec p = phi * particle[i].pbest_pos + (1 - phi) * pso::gbest_position;

		double k = random_double(0, 1, g);
    	double u = random_double(0, 1, g);

		if(k >= 0.5) {
			particle[i].current_position[0] = p[0] + beta * std::fabs(mbest[0] - particle[i].current_position[0]) * log(1.0/u);
			particle[i].current_position[1] = p[1] + beta * std::fabs(mbest[1] - particle[i].current_position[1]) * log(1.0/u);
			particle[i].current_position[2] = p[2] + beta * std::fabs(mbest[2] - particle[i].current_position[2]) * log(1.0/u);
    	}
    	else {
			particle[i].current_position[0] = p[0] - beta * std::fabs(mbest[0] - particle[i].current_position[0]) * log(1.0/u);
			particle[i].current_position[1] = p[1] - beta * std::fabs(mbest[1] - particle[i].current_position[1]) * log(1.0/u);
			particle[i].current_position[2] = p[2] - beta * std::fabs(mbest[2] - particle[i].current_position[2]) * log(1.0/u);
    	}
		
		//give a random position, if outside the search box
		if(particle[i].current_position[0] < corner1[0])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[1] < corner1[1])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[2] < corner1[2])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		
		if(particle[i].current_position[0] > corner2[0])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[1] > corner2[1])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[2] > corner2[2])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);

	}

    void pso::computeNewOrientation(int i, int max_iter, int current_iter, rng& generator)
	{
		double beta = (1 - 0.5) * (max_iter - current_iter) / static_cast<fl>(max_iter) + 0.5;

		int num = pso::number;
		qt mbest_tmp(0.0, 0.0, 0.0, 0.0);
		for(int j = 1; j <= num; ++j) {
			mbest_tmp += particle[j].pbest_orientation;
		}
   		qt mbest = mbest_tmp / static_cast<double>(num);

		double phi = random_double(0, 1, g);
		qt p =  phi * particle[i].pbest_orientation + (1 - phi) * pso::gbest_orientation;

		double k = random_double(0, 1, g);
    	double u = random_double(0, 1, g);

		qt tmp_orientation = beta * (mbest - particle[i].current_orientation) * log(1.0/u);
		tmp_orientation = qt(std::fabs(tmp_orientation.R_component_1()),
                     		 std::fabs(tmp_orientation.R_component_2()),
                     		 std::fabs(tmp_orientation.R_component_3()),
                     		 std::fabs(tmp_orientation.R_component_4()));

		vec tmp_rotation = quaternion_to_angle(tmp_orientation);

		std::cout << "更新前的方向" << particle[i].current_orientation << std::endl;

		if(k >= 0.5)
		{
        	quaternion_increment(p, tmp_rotation);
			particle[i].current_orientation = p;

    	}
		else
		{
        	tmp_rotation -= 2 * tmp_rotation;       // 取相反数
        	quaternion_increment(p, tmp_rotation);
			particle[i].current_orientation = p;
    	}
		std::cout << "更新后的方向" << particle[i].current_orientation << std::endl;
	}

	void pso::computeNewTorsion(int i, rng& generator, sz which, int max_iter, int current_iter)
	{
		double beta = (1 - 0.5) * (max_iter - current_iter) / static_cast<fl>(max_iter) + 0.5;

		int num = pso::number;
		double mbest_tmp = 0.0;
		for(int j = 0; j < num; ++j) 
		{
        	mbest_tmp += particle[j].pbest_torsion[which];
    	}
		double mbest = mbest_tmp / static_cast<double>(num);

		double phi = random_double(0, 1, generator);
		double p = phi * particle[i].pbest_torsion[which] + (1 - phi) * pso::gbest_torsion[which];

		double k = random_double(0, 1, generator);
		double u = random_double(0, 1, generator);
		for (int j = 0; j < which; j++) {
			std::cout << "第" << j + 1 << "个扭转：更新前" << particle[i].current_torsion[j] << std::endl;
			if(k >= 0.5) {
        		particle[i].current_torsion[j] = p + beta * std::fabs(mbest - particle[i].current_torsion[j]) * log(1.0/u);
    		}
			else {
        		particle[i].current_torsion[j] = p - beta * std::fabs(mbest - particle[i].current_torsion[j]) * log(1.0/u);
    		}
			std::cout << "第" << j + 1 << "个扭转：更新后" << particle[i].current_torsion[j] << std::endl;
		}

	}
	
	
	void pso::updatePersonalBest(int i,double e)
	{
		particle[i].pbest_fit = e;
	}

	
	void pso::updateGlobalBest(int i)
	{
		pso::gbest_fit = particle[i].pbest_fit;
	}
	
	void pso::updateGlobalBest_1(fl en)
	{
		pso::gbest_fit = en;
	}

	
	double pso::getPersonalBest(int i)
	{
		return particle[i].pbest_fit;
	}

	void pso::updateBestPosition(int i,vec pos)
	{
		particle[i].pbest_pos = pos;
		
	}
	
	void pso::updateBestOrientation(int i, qt orientation)
	{
		particle[i].pbest_orientation = orientation;
	}
	
	void pso::updateBestTorsion(int i, fl torsion,sz which)
	{
		particle[i].pbest_torsion[which] = torsion;
	}

	
	vec pso::getCurrentPosition(int i)
	{
		return particle[i].current_position;
	}
	
	qt pso::getCurrentOrientation(int i)
	{
		return particle[i].current_orientation;
	}
	
	fl pso::getCurrentTorsion(int i,sz which)
	{
		return particle[i].current_torsion[which];
	}

