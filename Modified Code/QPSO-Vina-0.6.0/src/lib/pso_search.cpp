/*
        PSOVina version 2.0

        This file is revised from monte_carlo.cpp in AutoDock Vina.

        Authors: Giotto H. K. TAI  <giottotai@yahoo.com.hk>

                 Shirley W. I. SIU <shirleysiu@umac.mo>


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

#include "pso_search.h"
#include "coords.h"
#include "quasi_newton.h"
#include "pso_mutate.h"
#include "mutate.h"

#include <string>
#include <sstream>

#include <iostream>
#include <exception>
#include <vector>
#include <cmath>
#include <boost/program_options.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/thread/thread.hpp>
// #include "parallel_pso.h"  不能包含这一行
#include "parse_pdbqt.h"
#include "file.h"
#include "cache.h"
#include "non_cache.h"
#include "naive_non_cache.h"
#include "parse_error.h"
#include "everything.h"
#include "weighted_terms.h"
#include "current_weights.h"
#include "quasi_newton.h"
#include "tee.h"
#include "coords.h"

// void exportForMDSimulation(model m, model ref_m, output_type candidate, std::string& remark, int step) {
//     std::ostringstream oss;
// 	oss << "output_" << step << ".pdbqt";
// 	std::string output_tmp = oss.str();

//     m.set(candidate.c);
//     const fl lb = m.rmsd_lower_bound(ref_m);
//     const fl ub = m.rmsd_upper_bound(ref_m);

//     remark = vina_remark(candidate.e, lb, ub);

//     // 输出为pdbqt文件
//     ofile f(make_path(output_tmp));
//     m.write_model(f, 1);
// }

// output_type importMDSimulationResults(model& m) {
//     // wt
//     fl weight_gauss1      = -0.035579;
//     fl weight_gauss2      = -0.005156;
//     fl weight_repulsion   =  0.840245;
//     fl weight_hydrophobic = -0.035069;
//     fl weight_hydrogen    = -0.587439;
//     fl weight_rot         =  0.05846;

//     flv weights;
//     weights.push_back(weight_gauss1);
//     weights.push_back(weight_gauss2);
//     weights.push_back(weight_repulsion);
//     weights.push_back(weight_hydrophobic);
//     weights.push_back(weight_hydrogen);
//     weights.push_back(5 * weight_rot / 0.1 - 1); 

//     everything t;
//     weighted_terms wt(&t, weights);

//     // prec
//     precalculate prec(wt);

//     // nnc
//     naive_non_cache nnc(&prec);

//     // authentic_v
//     const vec authentic_v(1000, 1000, 1000);

//     // c
//     conf c = m.get_initial_conf();

//     // intramolecular_energy
//     fl intramolecular_energy = m.eval_intramolecular(prec, authentic_v, c);
//     std::string ligand_name;
//     std::cout << "请输入ligand文件名：" << std::endl;
//     std::cin >> ligand_name;

//     std::cout << ligand_name << std::endl;
//     // 读取文件产生新的model
//     m = parse_bundle(std::vector<std::string>(1, ligand_name));

//     // 通过c和e产生新的candidate
//     fl e = m.eval_adjusted(wt, prec, nnc, authentic_v, c, intramolecular_energy);
//     output_type candidate(c, e);

//     return candidate;
// }

output_type pso_search::operator()(model& m, model ref_m, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator,int num_birds, double w, double c1, double c2) const {
	output_container tmp;
	this->operator()(m, ref_m, tmp, p, ig, p_widened, ig_widened, corner1, corner2, increment_me, generator, num_birds, w, c1, c2); // call the version that produces the whole container
	VINA_CHECK(!tmp.empty());
	return tmp.front();
}

bool metropolis_accept(fl old_f, fl new_f, fl temperature, rng& generator) {
	if(new_f < old_f) return true;
	const fl acceptance_probability = std::exp((old_f - new_f) / temperature);
	return random_fl(0, 1, generator) < acceptance_probability;
}

// out is sorted
void pso_search::operator()(model& m, model ref_m, output_container& out, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator,int num_birds,double w,double c1,double c2) const {
	vec authentic_v(1000, 1000, 1000); // FIXME? this is here to avoid max_fl/max_fl
	conf_size s = m.get_size();
	change g(s);
	output_type tmp(s, 0);
	tmp.c.randomize(corner1, corner2, generator);  //first randomize
	fl best_e = max_fl;
	quasi_newton quasi_newton_par; quasi_newton_par.max_steps = ssd_par.evals;
	output_type tmp_rough = tmp;
	pso particle(num_birds, w, c1, c2,corner1,corner2,generator,tmp.c);
	double* PersonalBest = new double[100];
	for(int cou = 0; cou < 100; cou++)
		PersonalBest[cou] = 0;
	double energy = 0;
	int count = 0;

    // std::string outputInfo = "outputInfo.csv";

    // ofile file(make_path(outputInfo));
    // // std::ofstream file(outputInfo);
    // file << "hello\n";
    // file.flush();
    // std::string remark;

	VINA_U_FOR(step, num_steps) {
		if(increment_me)
			++(*increment_me);
		output_type candidate = tmp_rough;
		output_type candidate_1 = tmp;

		pso_mutate_conf(candidate, candidate_1, m, ref_m, mutation_amplitude, generator, &particle, PersonalBest, p, ig, g, hunt_cap, quasi_newton_par, step, num_steps); //for each particle loop

		// exportForMDSimulation(m, ref_m, candidate, remark, step);
        // file << step << remark << '\n';
        // file.flush();
        // std::cout << "写入：" << "第" << step << "次迭代的信息" << std::endl;

		// candidate = importMDSimulationResults(m);

		tmp_rough = candidate;
		if(step == 0 || metropolis_accept(tmp.e, candidate_1.e, temperature, generator))
		{
			tmp = candidate_1;
			m.set(tmp.c); // FIXME? useless?
			if(tmp.e < best_e || out.size() < num_saved_mins) {
				quasi_newton_par(m, p, ig, tmp, g, authentic_v);
				m.set(tmp.c); // FIXME? useless?
				tmp.coords = m.get_heavy_atom_movable_coords();
				add_to_output_container(out, tmp, min_rmsd, num_saved_mins); // 20 - max size
				if(tmp.e < best_e)
						best_e = tmp.e;
			}
		}

		// FIXME only for very promising ones

		/***Criteria defined by PSOVina***/
		if(std::abs(pso::gbest_fit - energy) < 0.0001)
		{
			count += 1;
			if(count > 350)
			{
				step = num_steps; //break the loop
				count = 0;
			}
		}else{
			energy = pso::gbest_fit;
			count =0;
		}
	
	}

	
	VINA_CHECK(!out.empty());
	VINA_CHECK(out.front().e <= out.back().e); // make sure the sorting worked in the correct order
}

//
