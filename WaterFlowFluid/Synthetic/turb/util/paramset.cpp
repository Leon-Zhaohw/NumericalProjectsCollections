/******************************************************************************
 * Parameter storage class
 *
 * Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys.
 *
 *****************************************************************************/

// paramset.cpp*
#define CORE_SOURCE
#include "globals.h"
#include "paramset.h"

// ParamSet Macros
#define ADD_PARAM_TYPE(T, vec) \
	(vec).push_back(new ParamSetItem<T>(name, (const T *)data, nItems))

#define LOOKUP_PTR(vec) \
	for (size_t i = 0; i < (vec).size(); ++i) \
		if ((vec)[i]->name == name) { \
			*nItems = (vec)[i]->nItems; \
			(vec)[i]->lookedUp = true; \
			return (vec)[i]->data; \
		} \
	return NULL
#define LOOKUP_ONE(vec) \
	for (size_t i = 0; i < (vec).size(); ++i) { \
		if ((vec)[i]->name == name && \
			(vec)[i]->nItems == 1) { \
			(vec)[i]->lookedUp = true; \
			return *((vec)[i]->data); \
}		} \
	return d

// ParamSet Methods
ParamSet::ParamSet(const ParamSet &p2) {
	*this = p2;
}

void ParamSet::AddSet(const ParamSet &p2) {
	unsigned i;	
	for (i = 0; i < p2.ints.size(); ++i)
		AddInt(p2.ints[i]->name, p2.ints[i]->data);
	for (i = 0; i < p2.bools.size(); ++i)
		AddBool(p2.bools[i]->name, p2.bools[i]->data);
	for (i = 0; i < p2.floats.size(); ++i)
		AddFloat(p2.floats[i]->name, p2.floats[i]->data);
	for (i = 0; i < p2.vectors.size(); ++i)
		AddVector(p2.vectors[i]->name, p2.vectors[i]->data);
	for (i = 0; i < p2.strings.size(); ++i)
		AddString(p2.strings[i]->name, p2.strings[i]->data);
}

ParamSet &ParamSet::operator=(const ParamSet &p2) {
	if (&p2 != this) {
		Clear(); unsigned i;
		for (i = 0; i < p2.ints.size(); ++i)
			ints.push_back(p2.ints[i]->Clone());
		for (i = 0; i < p2.bools.size(); ++i)
			bools.push_back(p2.bools[i]->Clone());
		for (i = 0; i < p2.floats.size(); ++i)
			floats.push_back(p2.floats[i]->Clone());
		for (i = 0; i < p2.vectors.size(); ++i)
			vectors.push_back(p2.vectors[i]->Clone());
		for (i = 0; i < p2.strings.size(); ++i)
			strings.push_back(p2.strings[i]->Clone());
	}
	return *this;
}

void ParamSet::override(const ParamSet &p2) {
	// note - this only overrides parameters that already
	// exist in the paramset!

	if (&p2 != this) {
		for (size_t i = 0; i < p2.ints.size(); ++i) {
			for (size_t j = 0; j < this->ints.size(); ++j) {
				if (p2.ints[i]->name.compare( this->ints[j]->name )==0) {
					fprintf(stderr,"Overriding int '%s' \n", this->ints[j]->name.c_str() );
					this->ints[j] = p2.ints[i]->Clone();
				}
			}
		}
		for (size_t i = 0; i < p2.bools.size(); ++i) {
			for (size_t j = 0; j < this->bools.size(); ++j) {
				if (p2.bools[i]->name.compare( this->bools[j]->name )==0) {
					fprintf(stderr,"Overriding bool '%s' \n", this->bools[j]->name.c_str() );
					this->bools[j] = p2.bools[i]->Clone();
				}
			}
		}
		for (size_t i = 0; i < p2.floats.size(); ++i) {
			for (size_t j = 0; j < this->floats.size(); ++j) {
				if (p2.floats[i]->name.compare( this->floats[j]->name )==0) {
					fprintf(stderr,"Overriding float '%s' \n", this->floats[j]->name.c_str() );
					this->floats[j] = p2.floats[i]->Clone();
				}
			}
		}
		for (size_t i = 0; i < p2.vectors.size(); ++i) {
			for (size_t j = 0; j < this->vectors.size(); ++j) {
				if (p2.vectors[i]->name.compare( this->vectors[j]->name )==0) {
					fprintf(stderr,"Overriding vec '%s' \n", this->vectors[j]->name.c_str() );
					this->vectors[j] = p2.vectors[i]->Clone();
				}
			}
		}
		for (size_t i = 0; i < p2.strings.size(); ++i) {
			for (size_t j = 0; j < this->strings.size(); ++j) {
				if (p2.strings[i]->name.compare( this->strings[j]->name )==0) {
					fprintf(stderr,"Overriding string '%s' \n", this->strings[j]->name.c_str() );
					this->strings[j] = p2.strings[i]->Clone();
				}
			}
		}
	}
	return;
}

void ParamSet::AddFloat(const string &name,
			            const float *data,
						int nItems) {
	EraseFloat(name);
	floats.push_back(new ParamSetItem<float>(name,
											 data,
											 nItems));
}
void ParamSet::AddInt(const string &name, const int *data, int nItems) {
	EraseInt(name);
	ADD_PARAM_TYPE(int, ints);
}
void ParamSet::AddBool(const string &name, const bool *data, int nItems) {
	EraseInt(name);
	ADD_PARAM_TYPE(bool, bools);
}

void ParamSet::AddVector(const string &name, const DDF::Vec3 *data, int nItems) {
	EraseVector(name);
	ADD_PARAM_TYPE(DDF::Vec3, vectors);
}

void ParamSet::AddString(const string &name, const string *data, int nItems) {
	EraseString(name);
	ADD_PARAM_TYPE(string, strings);
}

bool ParamSet::EraseInt(const string &n) {
	for (size_t i = 0; i < ints.size(); ++i)
		if (ints[i]->name == n) {
			delete ints[i];
			ints.erase(ints.begin() + i);
			return true;
		}
	return false;
}
bool ParamSet::EraseBool(const string &n) {
	for (size_t i = 0; i < bools.size(); ++i)
		if (bools[i]->name == n) {
			delete bools[i];
			bools.erase(bools.begin() + i);
			return true;
		}
	return false;
}

bool ParamSet::EraseFloat(const string &n) {
	for (size_t i = 0; i < floats.size(); ++i)
		if (floats[i]->name == n) {
			delete floats[i];
			floats.erase(floats.begin() + i);
			return true;
		}
	return false;
}

bool ParamSet::EraseVector(const string &n) {
	for (size_t i = 0; i < vectors.size(); ++i)
		if (vectors[i]->name == n) {
			delete vectors[i];
			vectors.erase(vectors.begin() + i);
			return true;
		}
	return false;
}

bool ParamSet::EraseString(const string &n) {
	for (size_t i = 0; i < strings.size(); ++i)
		if (strings[i]->name == n) {
			delete strings[i];
			strings.erase(strings.begin() + i);
			return true;
		}
	return false;
}
float ParamSet::FindOneFloat(const string &name,
                             float d) const {
	for (size_t i = 0; i < floats.size(); ++i)
		if (floats[i]->name == name &&
			floats[i]->nItems == 1) {
			floats[i]->lookedUp = true;
			return *(floats[i]->data);
		}
	return d;
}
const float *ParamSet::FindFloat(const string &name,
		int *nItems) const {
	for (size_t i = 0; i < floats.size(); ++i)
		if (floats[i]->name == name) {
			*nItems = floats[i]->nItems;
			floats[i]->lookedUp = true;
			return floats[i]->data;
		}
	return NULL;
}
const int *ParamSet::FindInt(const string &name, int *nItems) const {
	LOOKUP_PTR(ints);
}
const bool *ParamSet::FindBool(const string &name, int *nItems) const {
	LOOKUP_PTR(bools);
}
int ParamSet::FindOneInt(const string &name, int d) const {
	LOOKUP_ONE(ints);
}
bool ParamSet::FindOneBool(const string &name, bool d) const {
	LOOKUP_ONE(bools);
}
const DDF::Vec3 *ParamSet::FindVector(const string &name, int *nItems) const {
	LOOKUP_PTR(vectors);
}
DDF::Vec3 ParamSet::FindOneVector(const string &name, const DDF::Vec3 &d) const {
	LOOKUP_ONE(vectors);
}
const string *ParamSet::FindString(const string &name, int *nItems) const {
	LOOKUP_PTR(strings);
}
string ParamSet::FindOneString(const string &name, const string &d) const {
	LOOKUP_ONE(strings);
}

void ParamSet::ReportUnused() const {
  bool haveUnused = false;

#define CHECK_UNUSED(v, TYPE) \
	for (i = 0; i < (v).size(); ++i) \
		if (!(v)[i]->lookedUp) { haveUnused = true; \
			debMsg("PARAMSET-WARN","Parameter \""<<  (v)[i]->name.c_str() <<"\" not used, type=" TYPE); \
		}
	size_t i;
	CHECK_UNUSED(ints,"int");       
	CHECK_UNUSED(bools, "bool");
	CHECK_UNUSED(floats, "float");   
	CHECK_UNUSED(vectors, "vec3"); 
	CHECK_UNUSED(strings, "string");

	if(haveUnused) {
	  errFatal("ParamSet::ReportUnused","Stopping...", SIMWORLD_ERRPARSE );
	}
}
void ParamSet::Clear() {
#define DEL_PARAMS(name) \
	for (size_t i = 0; i < (name).size(); ++i) \
		delete (name)[i]; \
	(name).erase((name).begin(), (name).end())
	DEL_PARAMS(ints);    
	DEL_PARAMS(bools);
	DEL_PARAMS(floats);  
	DEL_PARAMS(vectors); 
	DEL_PARAMS(strings);
#undef DEL_PARAMS
}
string ParamSet::ToString() const {
	string ret;
	size_t i;
	int j;
	string typeString;
	const int bufLen = 48*1024*1024;
	static char *buf = new char[bufLen];
	char *bufEnd = buf + bufLen;

	for (i = 0; i < ints.size(); ++i) {
		char *bufp = buf;
		*bufp = '\0';
		ParamSetItem<int> *item = ints[i];
		typeString = "integer ";
		// Print _ParamSetItem_ declaration, determine how many to print
		int nPrint = item->nItems;
		ret += string("\"");
		ret += typeString;
		ret += item->name;
		ret += string("\"");
		ret += string(" [");
		for (j = 0; j < nPrint; ++j)
			bufp += snprintf(bufp, bufEnd - bufp, "%d ", item->data[j]);
		ret += buf;
		ret += string("] ");
	}

	for (i = 0; i < bools.size(); ++i) {
		char *bufp = buf;
		*bufp = '\0';
		ParamSetItem<bool> *item = bools[i];
		typeString = "bool ";
		// Print _ParamSetItem_ declaration, determine how many to print
		int nPrint = item->nItems;
		ret += string("\"");
		ret += typeString;
		ret += item->name;
		ret += string("\"");
		ret += string(" [");
		for (j = 0; j < nPrint; ++j)
			bufp += snprintf(bufp, bufEnd - bufp, "\"%s\" ", item->data[j] ? "true" : "false");
		ret += buf;
		ret += string("] ");
	}

	for (i = 0; i < floats.size(); ++i) {
		char *bufp = buf;
		*bufp = '\0';
		ParamSetItem<float> *item = floats[i];
		typeString = "float ";
		// Print _ParamSetItem_ declaration, determine how many to print
		int nPrint = item->nItems;
		ret += string("\"");
		ret += typeString;
		ret += item->name;
		ret += string("\"");
		ret += string(" [");
		for (j = 0; j < nPrint; ++j)
			bufp += snprintf(bufp, bufEnd - bufp, "%.8g ", item->data[j]);
		ret += buf;
		ret += string("] ");
	}

	for (i = 0; i < vectors.size(); ++i) {
		char *bufp = buf;
		*bufp = '\0';
		ParamSetItem<DDF::Vec3> *item = vectors[i];
		typeString = "vector ";
		// Print _ParamSetItem_ declaration, determine how many to print
		int nPrint = item->nItems;
		ret += string("\"");
		ret += typeString;
		ret += item->name;
		ret += string("\"");
		ret += string(" [");
		for (j = 0; j < nPrint; ++j)
			bufp += snprintf(bufp, bufEnd - bufp, "%.8g %.8g %.8g ", item->data[j].x,
				item->data[j].y, item->data[j].z);
		ret += buf;
		ret += string("] ");
	}

	for (i = 0; i < strings.size(); ++i) {
		char *bufp = buf;
		*bufp = '\0';
		ParamSetItem<string> *item = strings[i];
		typeString = "string ";
		// Print _ParamSetItem_ declaration, determine how many to print
		int nPrint = item->nItems;
		ret += string("\"");
		ret += typeString;
		ret += item->name;
		ret += string("\"");
		ret += string(" [");
		for (j = 0; j < nPrint; ++j)
			bufp += snprintf(bufp, bufEnd - bufp, "\"%s\" ", item->data[j].c_str());
		ret += buf;
		ret += string("] ");
	}
	return ret;
}



