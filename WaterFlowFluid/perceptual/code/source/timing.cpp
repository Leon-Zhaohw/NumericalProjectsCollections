/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Plugin timing
 *
 ******************************************************************************/

#include "timing.h"
#include <fstream>
#include <iomanip>

using namespace std;
namespace Manta {

TimingData::TimingData() : updated(false), num(0) {
}

void TimingData::start(FluidSolver* parent, const string& name) {
	mLastPlugin = name;
	mPluginTimer.get();
}

void TimingData::stop(FluidSolver* parent, const string& name) {
	if (mLastPlugin == name && name != "FluidSolver::step") {
		updated = true;
		const string parentName = parent ? parent->getName() : "";
		MuTime diff = mPluginTimer.update();
		vector<TimingSet>& cur = mData[name];
		for (vector<TimingSet>::iterator it = cur.begin(); it != cur.end(); it++) {
			if (it->solver == parentName) {
				it->cur += diff;
				it->updated = true;
				return;
			}
		}
		TimingSet s;
		s.solver = parentName;
		s.cur = diff;
		s.updated = true;
		cur.push_back(s);
	}
}

void TimingData::step() {
	if (updated)
		num++;
	std::map<std::string, std::vector<TimingSet> >::iterator it;
	for (it = mData.begin(); it != mData.end(); it++) {
		for (vector<TimingSet>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
			if (it2->updated) {
				it2->total += it2->cur;
				it2->num++;
			}
			it2->cur.clear();
			it2->updated = false;
		}
	}
	updated = false;
}
 
void TimingData::print() {
	MuTime total;
	total.clear();
	std::map<std::string, std::vector<TimingSet> >::iterator it;
	for (it = mData.begin(); it != mData.end(); it++)
		for (vector<TimingSet>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++)
			total += it2->cur;

	std::cout << std::endl << "-- STEP " << std::setw(3) << num << " ----------------------------" << std::endl;
	for (it = mData.begin(); it != mData.end(); it++) {
		for (vector<TimingSet>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
			if (!it2->updated) continue;
			string name = it->first;
			if (it->second.size() > 1 && !it2->solver.empty())
				name += "[" + it2->solver + "]";
			std::cout << "["<< std::fixed << std::setw(5) << std::setprecision(2) <<
				(100.0*((Real)it2->cur.time / (Real)total.time)) << "%] " <<
				name << " (" << std::setprecision(5) << it2->cur << ")" << std::endl;
		}
	}
	step();

	std::cout << std::string(40, '-') << std::endl <<
		"Total: " << std::fixed << std::setprecision(5) << total << std::endl << std::endl;
}

void TimingData::saveMean(const string& filename) {
	ofstream ofs(filename.c_str());
	if (!ofs.good()) errMsg("can't open " + filename + " as timing log");
	ofs << "Mean timings of " << num << " steps :" <<endl <<endl;
	MuTime total;
	total.clear();
	std::map< std::string, std::vector<TimingSet> >::iterator it;
	for (it = mData.begin(); it != mData.end(); it++)
		for (vector<TimingSet>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
			total += it2->total;
			string name = it->first;
			if (it->second.size() > 1)
				name += "[" + it2->solver + "]";

			ofs << name << " " << (it2->total / it2->num) << " (" << it2->num << " calls)" << endl;
		}
	 
	ofs << endl << "Total : " << total.toSecond() << "s (mean " << total.toSecond()/num << "s)" << endl;
	ofs.close();
}
 
}

