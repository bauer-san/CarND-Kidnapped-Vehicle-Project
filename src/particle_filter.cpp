/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles=10; //1000;
	particles.resize(num_particles);

	default_random_engine gen;

	// from L14.5
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	for (int i=0;i<num_particles;i++){
		particles[i].x=dist_x(gen);
		particles[i].y=dist_y(gen);
		particles[i].theta=dist_theta(gen);
		particles[i].weight=1.;
	}

	is_initialized=true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle USING A MOTION MODEL and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;

	normal_distribution<double> dist_x(0., std_pos[0]);
	normal_distribution<double> dist_y(0., std_pos[1]);
	normal_distribution<double> dist_theta(0., std_pos[2]);

	for (int i=0;i<num_particles;i++){
		if (abs(yaw_rate)>0.001) {
			particles[i].x+=velocity/yaw_rate*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta)) + dist_x(gen);
			particles[i].y+=velocity/yaw_rate*(cos(particles[i].theta)-cos(particles[i].theta+yaw_rate*delta_t)) + dist_y(gen);
			particles[i].theta+=+yaw_rate*delta_t + dist_theta(gen);
		} else { // => CONSTANT yaw_angle
			particles[i].x+=velocity*delta_t*cos(particles[i].theta) + dist_x(gen);
			particles[i].y+=velocity*delta_t*sin(particles[i].theta) + dist_y(gen);
		} /*if*/
	}/*for*/
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	struct TransformedObs {
		double x;
		double y;
	} TObs;
	struct NN {
		int ID;
		double dist_m;
	} NearestNeighbor;

	double PtoL_dist;

	for (int i=0;i<num_particles;i++){
		for (int j=0;j<observations.size();j++){
			// Transform observations from VEHICLE to MAP coordinates
			TObs.x=observations[j].x*cos(particles[i].theta) - observations[j].y*sin(particles[i].theta) + particles[i].x;
			TObs.y=observations[j].x*sin(particles[i].theta) + observations[j].y*cos(particles[i].theta) + particles[i].y,

			NearestNeighbor.ID=0;
			NearestNeighbor.dist_m=9999.; // Initialize to a huge value
			for (int k=0;k<map_landmarks.landmark_list.size();k++){
				//Find distance from particle to Landmark
				PtoL_dist = dist(TObs.x, TObs.y, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f);

				//If a new NN was found update distance and ID
				if (PtoL_dist<NearestNeighbor.dist_m){
					NearestNeighbor.dist_m=PtoL_dist;
					NearestNeighbor.ID=map_landmarks.landmark_list[k].id_i;
				}/*if*/
			}/*for k*/
		}/*for j*/
	}/*for i*/
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
