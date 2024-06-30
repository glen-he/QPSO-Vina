/*
        PSOVina version 1.0                     Date: 26/11/2014

        This file is revised from mutate.h in AutoDock Vina.

        Authors: Marcus C. K. Ng  <marcus.ckng@gmail.com>

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

#ifndef VINA_PSO_MUTATE_H
#define VINA_PSO_MUTATE_H

#include "pso.h"
#include "model.h"
#include "quasi_newton.h"


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

// does not set model
void pso_mutate_conf(output_type& c, output_type& c_1, model& m, model ref_m, fl amplitude, rng& generator, pso*, double*, const precalculate&, const igrid&, change&, const vec&, quasi_newton&, int, int);
void pso_mutate_conf(output_type& c, model& m, model ref_m, fl amplitude, rng& generator,pso*,const precalculate&,const igrid&,change&,const vec&,quasi_newton&,int, int);
#endif
