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

/******************************************************/
#include "pso_search.h"
#include "coords.h"
#include "quasi_newton.h"
#include "pso_mutate.h"
#include "mutate.h"

#include <string>
#include <sstream>
#include <iomanip>

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
/******************************************************/

/******************************************************/
std::string num_to_string(int num) {
    std::string result;
    bool isNegative = false;
    if (num < 0) {
        isNegative = true;
        num = -num;
    }
    do {
        result.insert(result.begin(), num % 10 + '0');
        num /= 10;
    } while (num != 0);
    if (isNegative) {
        result.insert(result.begin(), '-');
    }
    return result;
}

void exportForMDSimulation(model m, model ref_m, output_type candidate, std::string& remark, int step, int number) {
    std::ostringstream oss;
	oss << "dk_s" << std::setw(2) << std::setfill('0') << step + 1 << "p" << std::setw(2) << std::setfill('0') << number + 1<< ".pdbqt";
	std::string output_tmp = oss.str();

    m.set(candidate.c);
    const fl lb = m.rmsd_lower_bound(ref_m);
    const fl ub = m.rmsd_upper_bound(ref_m);

    remark = vina_remark(candidate.e, lb, ub);

    // 输出为pdbqt文件
    ofile f(make_path(output_tmp));
    m.write_model(f, 1);
}

output_type importMDSimulationResults(model& m, int step, int number) {
    // wt
    fl weight_gauss1      = -0.035579;
    fl weight_gauss2      = -0.005156;
    fl weight_repulsion   =  0.840245;
    fl weight_hydrophobic = -0.035069;
    fl weight_hydrogen    = -0.587439;
    fl weight_rot         =  0.05846;

    flv weights;
    weights.push_back(weight_gauss1);
    weights.push_back(weight_gauss2);
    weights.push_back(weight_repulsion);
    weights.push_back(weight_hydrophobic);
    weights.push_back(weight_hydrogen);
    weights.push_back(5 * weight_rot / 0.1 - 1); 

    everything t;
    weighted_terms wt(&t, weights);

    // prec
    precalculate prec(wt);

    // nnc
    naive_non_cache nnc(&prec);

    // authentic_v
    const vec authentic_v(1000, 1000, 1000);

    // c
    conf c = m.get_initial_conf();

    // intramolecular_energy
    fl intramolecular_energy = m.eval_intramolecular(prec, authentic_v, c);
    std::string ligand_name;
	std::cout << "输入第" << std::setw(2) << std::setfill('0') << step + 1 << "步，第" << std::setw(2) << std::setfill('0') << number + 1 << "个粒子，经过MD模拟后的配体" << std::endl;
    std::cin >> ligand_name;

    std::cout << "输入的信息为：" << ligand_name << std::endl;
    // 读取文件产生新的model
    m = parse_bundle(std::vector<std::string>(1, ligand_name));

    // 通过c和e产生新的candidate
    fl e = m.eval_adjusted(wt, prec, nnc, authentic_v, c, intramolecular_energy);
    output_type candidate(c, e);

    return candidate;
}



/******************************************************/

sz pso_count_mutable_entities(const conf& c) {
	sz counter = 0;
	VINA_FOR_IN(i, c.ligands)
		counter += 2 + c.ligands[i].torsions.size();
	VINA_FOR_IN(i, c.flex)
		counter += c.flex[i].torsions.size();
	return counter;
}

void pso_mutate_conf(output_type& candidate, output_type& candidate_1, model& m, model ref_m, fl amplitude, rng& generator, pso* particle, double* PersonalBest, const precalculate& p ,const igrid& ig,change& g,const vec& v,quasi_newton& quasi_newton_par,int step, int num_steps) { // ONE OF: 2A for position, similar amp for orientation, randomize torsion
	sz mutable_entities_num = pso_count_mutable_entities(candidate.c);
	if(mutable_entities_num == 0) return;
	// int which_int = random_int(0, int(mutable_entities_num - 1), generator);   // GLEN：不再使用随机的数字
	int which_int = mutable_entities_num;  // GLEN：改为使用所有的可变实体
	sz which = sz(which_int);

	std::cout << "which的值" << which << std::endl;
	
	srand((unsigned)time(NULL));
	double r_cm;
	r_cm = (double) rand() / (RAND_MAX + 1.0);
	
	double w_cm;
	w_cm = (double) rand() / (RAND_MAX + 1.0);

	int y;

	/* 输出Docking之后信息 */
	std::string outputInfo_Docking = "OutInfo_Docking_s" + num_to_string(step + 1) + ".csv";

	ofile file_dk(make_path(outputInfo_Docking));
	file_dk << "记录Docking后每个粒子信息\n";
	file_dk.flush();
	std::string remark;

	/* 输出MD之后信息 */
	std::string outputInfo_MD = "OutInfo_MD_s" + num_to_string(step + 1) + ".csv";

	ofile file_md(make_path(outputInfo_MD));
	file_md << "记录MD后的能量值\n";
	file_md.flush();

	VINA_FOR_IN (i, candidate.c.ligands) {
		model tmp_m = m;
		const vec authentic_v(1000, 1000, 1000);
		fl gr = m.gyration_radius(i);  // gr是配体节点原点与附着原子距离平方的平均值的平方根

		/* 分为两种情况，如果gr > epsion_fl，进行位置、方向、扭转的更新，反之只进行位置和扭转的更新 */

		// 情况1：gr > epsion_fl
		if (gr > epsilon_fl) {  // epsilon_fl是浮点数表示的精度限制

			for (y = 0; y < particle->number; y++) {    // 遍历每个粒子

				candidate.c.ligands[i].rigid.position = particle->getCurrentPosition(y);    // 获取当前粒子的位置
				candidate.c.ligands[i].rigid.orientation = particle->getCurrentOrientation(y);    // 获取当前粒子的方向
				for(int z=0; z<candidate.c.ligands[i].torsions.size(); z++) {
					candidate.c.ligands[i].torsions[z] = particle->getCurrentTorsion(y,z);    // 获取当前粒子每一个可旋转角的扭转
				}

				/* 输出粒子进行MD模拟，然后输入回来 */
				exportForMDSimulation(m, ref_m, candidate, remark, step, y);    // 输出粒子用于MD模拟
				file_dk << "第" << std::setw(2) << std::setfill('0') << step + 1 << "次迭代的第" << std::setw(2) << std::setfill('0') << y + 1 << "个粒子的信息：	" << remark;    // 存储当前粒子的reamrk信息到OutputInfo文件中
				file_dk.flush();

				candidate = importMDSimulationResults(m, step, y);    // MD模拟后到粒子输入回来
				const fl lb_tmp = m.rmsd_lower_bound(ref_m);
    			const fl ub_tmp = m.rmsd_upper_bound(ref_m);
				file_md << "第" << std::setw(2) << std::setfill('0') << step + 1 << "次迭代的第" << std::setw(2) << std::setfill('0') << y + 1 << "个粒子MD后的能量值：	" << candidate.e << "	" << lb_tmp << "	" << ub_tmp << "\n";    // 存储当前粒子的reamrk信息到OutputInfo文件中
				file_md.flush();

				/* 如果MD模拟后到粒子能量小于当前的PBest */ 
				if (candidate.e < particle->getPersonalBest(y)) {

					PersonalBest[y] = candidate.e;    // 更新粒子能量的PBest

					particle->updateBestPosition(y, candidate.c.ligands[i].rigid.position);    // 更新粒子位置的PBest
					particle->updateBestOrientation(y, candidate.c.ligands[i].rigid.orientation);    // 更新粒子方向的PBest
					for (int z = 0; z < candidate.c.ligands[i].torsions.size(); z++) {
						particle->updateBestTorsion(y, candidate.c.ligands[i].torsions[z],z);    // 更新粒子扭转的PBest
					}

					/* 如果MD模拟后到粒子能量小于当前的GBest */ 
					if (candidate.e < pso::gbest_fit) {

						particle->updateGlobalBest_1(candidate.e);    // 更新粒子能量的GBest

						pso::gbest_position = candidate.c.ligands[i].rigid.position;    // 更新粒子位置的GBest
						pso::gbest_orientation = candidate.c.ligands[i].rigid.orientation;    // 更新粒子方向的GBest
						for(int z = 0; z < candidate.c.ligands[i].torsions.size(); z++) {
							pso::gbest_torsion[z] = candidate.c.ligands[i].torsions[z];    // 更新粒子扭转的GBest
						}
					}
				}

				/* QPSO位置、方向、扭转的更新 */
				r_cm = 1.07*(7.86*r_cm-23.31*r_cm*r_cm+28.75*r_cm*r_cm*r_cm-13.302875*r_cm*r_cm*r_cm*r_cm);
				w_cm = 1.07*(7.86*w_cm-23.31*w_cm*w_cm+28.75*w_cm*w_cm*w_cm-13.302875*w_cm*w_cm*w_cm*w_cm);

				particle->updateVelocity(generator,y,r_cm,w_cm);    // 位置的速度
				particle->updateVelocityO(generator,y,r_cm,w_cm);    // 方向的速度
				particle->updateVelocityT(generator, y, which - 2, r_cm, w_cm);    // 扭转的速度

				particle->computeNewPositions(y, num_steps + 1, step + 1, generator);    // 计算粒子新的位置
				particle->computeNewOrientation(y, num_steps + 1, step + 1, generator);    // 计算粒子新的方向
				particle->computeNewTorsion(y, generator, which - 2, num_steps + 1, step + 1);    // 计算粒子新的扭转
			}

			/* 记录最优的粒子并赋值给candidate_1 */
			candidate_1.e = pso::gbest_fit;    // 将能量的GBest赋值给candidate_1

			candidate_1.c.ligands[i].rigid.position = pso::gbest_position;    // 将位置的GBest赋值给candidate_1
			candidate_1.c.ligands[i].rigid.orientation = pso::gbest_orientation;    // 将方向的GBest赋值给candidate_1
			for(int z = 0; z < candidate.c.ligands[i].torsions.size(); z++) {
				candidate_1.c.ligands[i].torsions[z] = pso::gbest_torsion[z];    // 将扭转的GBest赋值给candidate_1
			}
		}

		// 情况2：gr < epsion_fl
		else {

			for (y = 0; y < particle->number; y++) {    // 遍历每个粒子

				candidate.c.ligands[i].rigid.position = particle->getCurrentPosition(y);    // 获取当前粒子的位置
				candidate.c.ligands[i].rigid.orientation = particle->getCurrentOrientation(y);    // 获取当前粒子的方向
				for(int z=0; z<candidate.c.ligands[i].torsions.size(); z++) {
					candidate.c.ligands[i].torsions[z] = particle->getCurrentTorsion(y,z);    // 获取当前粒子每一个可旋转角的扭转
				}

				/* 输出粒子进行MD模拟，然后输入回来 */
				exportForMDSimulation(m, ref_m, candidate, remark, step, y);    // 输出粒子用于MD模拟
				file_dk << "第" << std::setw(2) << std::setfill('0') << step + 1 << "次迭代的第" << std::setw(2) << std::setfill('0') << y + 1 << "个粒子的信息：	" << remark;    // 存储当前粒子的reamrk信息到OutputInfo文件中
				file_dk.flush();

				candidate = importMDSimulationResults(m, step, y);    // MD模拟后到粒子输入回来
				const fl lb_tmp = m.rmsd_lower_bound(ref_m);
    			const fl ub_tmp = m.rmsd_upper_bound(ref_m);
				file_md << "第" << std::setw(2) << std::setfill('0') << step + 1 << "次迭代的第" << std::setw(2) << std::setfill('0') << y + 1 << "个粒子MD后的能量值：	" << candidate.e << "	" << lb_tmp << "	" << ub_tmp << "\n";    // 存储当前粒子的reamrk信息到OutputInfo文件中
				file_md.flush();

				/* 如果MD模拟后到粒子能量小于当前的PBest */ 
				if (candidate.e < particle->getPersonalBest(y)) {

					PersonalBest[y] = candidate.e;    // 更新粒子能量的PBest

					particle->updateBestPosition(y, candidate.c.ligands[i].rigid.position);    // 更新粒子位置的PBest
					particle->updateBestOrientation(y, candidate.c.ligands[i].rigid.orientation);    // 更新粒子方向的PBest
					for (int z = 0; z < candidate.c.ligands[i].torsions.size(); z++) {
						particle->updateBestTorsion(y, candidate.c.ligands[i].torsions[z], z);    // 更新粒子扭转的PBest
					}

					/* 如果MD模拟后到粒子能量小于当前的GBest */ 
					if (candidate.e < pso::gbest_fit) {

						particle->updateGlobalBest_1(candidate.e);    // 更新粒子能量的GBest

						pso::gbest_position = candidate.c.ligands[i].rigid.position;    // 更新粒子位置的GBest
						pso::gbest_orientation = candidate.c.ligands[i].rigid.orientation;    // 更新粒子方向的GBest
						for(int z = 0; z < candidate.c.ligands[i].torsions.size(); z++) {
							pso::gbest_torsion[z] = candidate.c.ligands[i].torsions[z];    // 更新粒子扭转的GBest
						}
					}
				}

				/* QPSO位置、方向、扭转的更新 */
				r_cm = 1.07*(7.86*r_cm-23.31*r_cm*r_cm+28.75*r_cm*r_cm*r_cm-13.302875*r_cm*r_cm*r_cm*r_cm);
				w_cm = 1.07*(7.86*w_cm-23.31*w_cm*w_cm+28.75*w_cm*w_cm*w_cm-13.302875*w_cm*w_cm*w_cm*w_cm);

				particle->updateVelocity(generator, y, r_cm, w_cm);    // 位置的速度
				// particle->updateVelocityO(generator,y,r_cm,w_cm);    // 方向的速度  情况二不需要改变
				particle->updateVelocityT(generator, y, which - 2, r_cm,w_cm);    // 扭转的速度

				particle->computeNewPositions(y, num_steps + 1, step + 1, generator);    // 计算粒子新的位置
				// particle->computeNewOrientation(y, num_steps, step, generator);    // 计算粒子新的方向  情况二不需要改变
				particle->computeNewTorsion(y, generator, which - 2, num_steps + 1, step + 1);    // 计算粒子新的扭转
			}

			/* 记录最优的粒子并赋值给candidate_1 */
			candidate_1.e = pso::gbest_fit;    // 将能量的GBest赋值给candidate_1

			candidate_1.c.ligands[i].rigid.position = pso::gbest_position;    // 将位置的GBest赋值给candidate_1
			candidate_1.c.ligands[i].rigid.orientation = pso::gbest_orientation;    // 将方向的GBest赋值给candidate_1
			for(int z = 0; z < candidate.c.ligands[i].torsions.size(); z++) {
				candidate_1.c.ligands[i].torsions[z] = pso::gbest_torsion[z];    // 将扭转的GBest赋值给candidate_1
			}
		}
		return;
	}
	file_dk.close();
	file_md.close();

	VINA_FOR_IN(i, candidate.c.flex) {
		if(which < candidate.c.flex[i].torsions.size()) { candidate.c.flex[i].torsions[which] = random_fl(-pi, pi, generator); return; }
		which -= candidate.c.flex[i].torsions.size();
	}
}
