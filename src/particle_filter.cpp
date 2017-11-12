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

	num_particles=2;
	particles.resize(num_particles);

	default_random_engine gen;

	normal_distribution<double> NormDist_x(x, std[0]);
	normal_distribution<double> NormDist_y(y, std[1]);
	normal_distribution<double> NormDist_theta(theta, std[2]);

	for (int i=0;i<num_particles;i++){
		Particle this_particle;

		this_particle.id=i;
		this_particle.x=NormDist_x(gen);
		this_particle.y=NormDist_y(gen);
		this_particle.theta=NormDist_theta(gen);
		this_particle.weight=1.;

		particles.push_back(this_particle);
		weights.push_back(1.);
	}

	is_initialized=true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle USING A MOTION MODEL and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;

	// Use zero mean and given std
	normal_distribution<double> NormDist_x(0., std_pos[0]);
	normal_distribution<double> NormDist_y(0., std_pos[1]);
	normal_distribution<double> NormDist_theta(0., std_pos[2]);

	Particle this_particle;

	for (int i=0;i<num_particles;i++){
		this_particle = particles[i];
		if (yaw_rate==0.) { // => CONSTANT yaw_angle
			this_particle.x+=velocity*delta_t*cos(this_particle.theta) + NormDist_x(gen);
			this_particle.y+=velocity*delta_t*sin(this_particle.theta) + NormDist_y(gen);
		} else {
			this_particle.x+=velocity/yaw_rate*(sin(this_particle.theta+yaw_rate*delta_t)-sin(this_particle.theta)) + NormDist_x(gen);
			this_particle.y+=velocity/yaw_rate*(cos(this_particle.theta)-cos(this_particle.theta+yaw_rate*delta_t)) + NormDist_y(gen);
			this_particle.theta+=+yaw_rate*delta_t + NormDist_theta(gen);
		} /*if*/

	particles[i]=this_particle;

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

	//For each particle:
	//		For Each Observation:
	//			Transform observations from VEHICLE to MAP coordinates -> Calculate TObs
	//			Associate the nearest landmark with each TObs 
	//			Update Weight -> Apply Multivariate-Gaussian Probability to calculate P_xy
	//			Combine probabilities (multiply)

	LandmarkObs TObs;

	int lm_id;
	double lm_x, lm_y, mu_x, mu_y;

	double PtoL_dist;

	struct NN {
		int ID;
		double dist_m;
	} NearestNeighbor;

	std::vector<int> associations;
	std::vector<double> sense_x;
	std::vector<double> sense_y;

	double P_xy, scaleP_xy = 1./(2.*M_PI*std_landmark[0]*std_landmark[1]);  //1.76839

	bool hit;

	for (int i=0;i<num_particles;i++){ // for each PARTICLE
		//particles[i].weight=1.;
		hit = false;
		for (int j=0;j<observations.size();j++){ // for each OBSERVATION
			// Transform the observation from VEHICLE to MAP coordinates
			TObs.x=observations[j].x*cos(particles[i].theta) - observations[j].y*sin(particles[i].theta) + particles[i].x;
			TObs.y=observations[j].x*sin(particles[i].theta) + observations[j].y*cos(particles[i].theta) + particles[i].y,

			// Associate the nearest landmark with each TObs
			NearestNeighbor.ID=0;
			NearestNeighbor.dist_m=99999.; // Initialize to a huge value
			for (int k=0;k<map_landmarks.landmark_list.size();k++){
				lm_id = map_landmarks.landmark_list[k].id_i;
				lm_x = map_landmarks.landmark_list[k].x_f;
				lm_y = map_landmarks.landmark_list[k].y_f;

				// Calculate the distance from Observation to Landmark 
				PtoL_dist=dist(TObs.x, TObs.y, lm_x, lm_y);

				// If new NearestNeighbor, update distance and ID
				if (PtoL_dist<NearestNeighbor.dist_m){
					NearestNeighbor.dist_m=PtoL_dist;
					NearestNeighbor.ID=lm_id;
//cout << i << "\t" << NearestNeighbor.ID << endl;
				}/*if*/
			}/*for k*/

			mu_x = map_landmarks.landmark_list[NearestNeighbor.ID-1].x_f;
			mu_y = map_landmarks.landmark_list[NearestNeighbor.ID-1].y_f;

			// Update Weight -> Apply Multivariate-Gaussian Probability to calculate P_xy
			P_xy = scaleP_xy * exp(-(((TObs.x-mu_x)*(TObs.x-mu_x))/(2.*std_landmark[0]*std_landmark[0]) + ((TObs.y-mu_y)*(TObs.y-mu_y))/(2.*std_landmark[1]*std_landmark[1])));
//cout << NearestNeighbor.ID << "\t" << TObs.x << "\t" << mu_x << "\t" << TObs.y << "\t" << mu_y << "\tP_xy: " << P_xy << endl;
//cout << P_xy << endl;
			if (P_xy !=0) {
				if (!hit) {
					particles[i].weight = P_xy;	// Combine probabilities
					hit=true;
				} else {
					particles[i].weight *= P_xy;	// Combine probabilities
				}/*if !hit*/
			}/*if P_xy !=0*/
//cout << mu_x << "\t" << mu_y << "\t" << particles[i].weight << endl;
cout << "Particle: " << i << "\tWeight: " << particles[i].weight <<  "\tObservation: " << j << "\tLandmark: " << NearestNeighbor.ID << endl;
			// Update the association, sense_x, and sense_y for this observation
			associations.push_back(NearestNeighbor.ID);
			sense_x.push_back(mu_x);
			sense_y.push_back(mu_y);

		}/*for j*/

	particles[i] = SetAssociations(particles[i], associations, sense_x, sense_y);
	weights[i]=particles[i].weight;

	}/*for i*/
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;
	discrete_distribution<int> distribution(weights.begin(), weights.end());

	vector<Particle> resample_particles;

	for (int i=0;i<num_particles;i++){
		resample_particles.push_back(particles[distribution(gen)]);
	}

	particles=resample_particles;

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
