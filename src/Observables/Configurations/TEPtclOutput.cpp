/*
 * TEPtclOutput.h
 *
 *  Created on: 19/jun/2015
 *      Author: Ferdinando (after TcLOutput.cpp, written by Flavio)
 */

#include <sstream>
#include "TEPtclOutput.h"

template<typename number>
TEPtclOutput<number>::TEPtclOutput() : Configuration<number>() {
	_print_labels = false;
	_core_radius = 1;
	_side_radius = 0.7;
	_front_radius = 0.82;
	
	_side_shift=0.3;
	_front_shift=0.15;
	
	_resolution = 20;
	_ref_particle_id = -1;
	_ref_strand_id = -1;
	_back_in_box = true;

	OX_LOG(Logger::LOG_WARNING,"Observable TEPtclOutput has not been properly tested since vmd takes forever to read tcl commands. Use TEPxyzOutput instead.");
}

template<typename number>
TEPtclOutput<number>::~TEPtclOutput() {

}

template<typename number>
void TEPtclOutput<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	Configuration<number>::get_settings(my_inp, sim_inp);
	int tmp;
	if(getInputBoolAsInt(&my_inp, "back_in_box", &tmp, 0) == KEY_FOUND) _back_in_box = (bool)tmp;
	if(getInputBoolAsInt(&my_inp, "print_labels", &tmp, 0) == KEY_FOUND) _print_labels = (bool)tmp;
	if(getInputInt(&my_inp, "ref_particle", &tmp, 0) == KEY_FOUND) _ref_particle_id = tmp;
	if(getInputInt(&my_inp, "ref_strand", &tmp, 0) == KEY_FOUND) _ref_strand_id = tmp;
	if(getInputInt(&my_inp, "resolution", &tmp, 0) == KEY_FOUND) _resolution = tmp;
}

template<typename number>
std::string TEPtclOutput<number>::_headers(llint step) {
	std::stringstream headers;
	
	headers << "color Display Background white" << endl;
	headers << "mol new" << endl;

	// we might want to be able to change these in the future...
	double _box_radius = 0.1;
	int _box_resolution = _resolution;
	
	number mybox = this->_config_info.box->box_sides()[0];
	headers << "graphics 0 cylinder {" << -mybox / 2. << " " << -mybox / 2. << " " << -mybox / 2. << "} {" <<  mybox / 2. << " " << -mybox / 2. << " " << -mybox / 2. << "} radius " << _box_radius << " resolution " << _box_resolution << " filled yes" << endl;
	headers << "graphics 0 cylinder {" << -mybox / 2. << " " << -mybox / 2. << " " <<  mybox / 2. << "} {" <<  mybox / 2. << " " << -mybox / 2. << " " <<  mybox / 2. << "} radius " << _box_radius << " resolution " << _box_resolution << " filled yes" << endl;
	headers << "graphics 0 cylinder {" << -mybox / 2. << " " <<  mybox / 2. << " " << -mybox / 2. << "} {" << -mybox / 2. << " " << -mybox / 2. << " " << -mybox / 2. << "} radius " << _box_radius << " resolution " << _box_resolution << " filled yes" << endl;
	headers << "graphics 0 cylinder {" << -mybox / 2. << " " <<  mybox / 2. << " " <<  mybox / 2. << "} {" << -mybox / 2. << " " << -mybox / 2. << " " <<  mybox / 2. << "} radius " << _box_radius << " resolution " << _box_resolution << " filled yes" << endl;
	headers << "graphics 0 cylinder {" <<  mybox / 2. << " " <<  mybox / 2. << " " <<  mybox / 2. << "} {" <<  mybox / 2. << " " << -mybox / 2. << " " <<  mybox / 2. << "} radius " << _box_radius << " resolution " << _box_resolution << " filled yes" << endl;
	headers << "graphics 0 cylinder {" <<  mybox / 2. << " " <<  mybox / 2. << " " << -mybox / 2. << "} {" <<  mybox / 2. << " " << -mybox / 2. << " " << -mybox / 2. << "} radius " << _box_radius << " resolution " << _box_resolution << " filled yes" << endl;
	headers << "graphics 0 cylinder {" << -mybox / 2. << " " <<  mybox / 2. << " " << -mybox / 2. << "} {" <<  mybox / 2. << " " <<  mybox / 2. << " " << -mybox / 2. << "} radius " << _box_radius << " resolution " << _box_resolution << " filled yes" << endl;
	headers << "graphics 0 cylinder {" << -mybox / 2. << " " <<  mybox / 2. << " " <<  mybox / 2. << "} {" <<  mybox / 2. << " " <<  mybox / 2. << " " <<  mybox / 2. << "} radius " << _box_radius << " resolution " << _box_resolution << " filled yes" << endl;
	headers << "graphics 0 cylinder {" << -mybox / 2. << " " << -mybox / 2. << " " << -mybox / 2. << "} {" << -mybox / 2. << " " << -mybox / 2. << " " <<  mybox / 2. << "} radius " << _box_radius << " resolution " << _box_resolution << " filled yes" << endl;
	headers << "graphics 0 cylinder {" <<  mybox / 2. << " " << -mybox / 2. << " " << -mybox / 2. << "} {" <<  mybox / 2. << " " << -mybox / 2. << " " <<  mybox / 2. << "} radius " << _box_radius << " resolution " << _box_resolution << " filled yes" << endl;
	headers << "graphics 0 cylinder {" <<  mybox / 2. << " " <<  mybox / 2. << " " << -mybox / 2. << "} {" <<  mybox / 2. << " " <<  mybox / 2. << " " <<  mybox / 2. << "} radius " << _box_radius << " resolution " << _box_resolution << " filled yes" << endl;
	headers << "graphics 0 cylinder {" << -mybox / 2. << " " <<  mybox / 2. << " " << -mybox / 2. << "} {" << -mybox / 2. << " " <<  mybox / 2. << " " <<  mybox / 2. << "} radius " << _box_radius << " resolution " << _box_resolution << " filled yes" << endl;

	return headers.str();
}

template<typename number>
std::string TEPtclOutput<number>::_particle(BaseParticle<number> *p) {
	std::stringstream res;
	
	TEPParticle<number> * me;
	me = reinterpret_cast<TEPParticle<number> *> (p);
	//next = reinterpret_cast<TEPParticle<number> *> (p->n5);
	
	LR_vector<number> zero (0., 0., 0.);
	if (_ref_particle_id >= 0 && this->_visible_particles.count(_ref_particle_id) == 1) zero = this->_config_info.particles[_ref_particle_id]->pos; 
	if (_ref_strand_id >= 0 && this->_strands_cdm.count(_ref_strand_id) == 1) zero = this->_strands_cdm[_ref_strand_id];
		
	// set the colour according to the strand id
	int colorid = p->strand_id;
	if (colorid >= 8) colorid ++;
	colorid = colorid % 33;
	
	res << "graphics 0 color " << colorid << endl;
	BaseBox<number> * mybox = this->_config_info.box;
	LR_vector<number> my_strand_cdm = this->_strands_cdm[me->strand_id];
	LR_vector<number> origin (0., 0., 0.);
	origin = zero;

	LR_vector<number> core = (me->pos - my_strand_cdm) + mybox->min_image(origin, my_strand_cdm);
	LR_vector<number> side = (me->pos - my_strand_cdm) + mybox->min_image(origin, my_strand_cdm) + _side_shift*me->orientationT.v2; 
	LR_vector<number> front = (me->pos - my_strand_cdm) + mybox->min_image(origin, my_strand_cdm) + _front_shift*me->orientationT.v1; 

	
	if (_print_labels && p->n5 != P_VIRTUAL && p->n3 == P_VIRTUAL) res << "graphics 0 text {" << core.x << " " << core.y << " " << core.z << "} \"" << me->strand_id << "\" size 1.5 " << endl;
	
	// core
	res << "graphics 0 sphere {" << core.x << " " << core.y << " " << core.z << "} radius " << _core_radius << " resolution " << _resolution << endl;
	// base
	res << "graphics 0 sphere {" << side.x << " " << side.y << " " << side.z << "} radius " << _side_radius << " resolution " << _resolution << endl;
	// front 
	res << "graphics 0 sphere {" << front.x << " " << front.y << " " << front.z << "} radius " << _front_radius << " resolution " << _resolution << endl;
	
	
	res << "graphics 0 color white";

	return res.str();
}


template<typename number>
std::string TEPtclOutput<number>::_configuration(llint step) {
	stringstream conf;
	//conf.precision(15);
	if (_back_in_box) this->_fill_strands_cdm ();
	//return Configuration<number>::_configuration(step);
	for(set<int>::iterator it = this->_visible_particles.begin(); it != this->_visible_particles.end(); it++) {
		if(it != this->_visible_particles.begin()) conf << endl;
		BaseParticle<number> *p = this->_config_info.particles[*it];
		conf << _particle(p);
	}
	return conf.str();
}

template class TEPtclOutput<float>;
template class TEPtclOutput<double>;

