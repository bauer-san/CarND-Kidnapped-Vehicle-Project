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

	num_particles=20;
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
	//
	// INPUTS:
	//		predicted => what we expected to see (i.e., Map Landmarks - in particle reference frame)
	//		observations => what we actually see


	// For each Observation
	//		For each in-range Landmark
	//			calculate distance
	//		Assign ID of Landmark closest to this observation to Observation ID

	double min_dist;

	for (int j=0;j<observations.size();j++){ // for each OBSERVATION

		min_dist=-1.;  // RESET min_dist

		for (int k=0;k<predicted.size();k++){ // for each PREDICTED / LANDMARK
			double this_dist;
		
			this_dist=dist(predicted[k].x, predicted[k].y, observations[j].x, observations[j].y);

			// FIND THE NEAREST NEIGHBOR
			if ((this_dist<min_dist) || (min_dist==-1.)) {  // if first or new minimum distance
				min_dist=this_dist;
				observations[j].id=predicted[k].id;			// HELLO NEIGHBOR
			}
		} // for k (LANDMARK)
//cout << "\tObservation: " << j << "\tNearest Neighbor: " << observations[j].id << endl;

	} // for j OBSERVATION

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
	//	Get a list of all landmarks in sensor_range
	//		For Each Observation:
	//			Transform observations from VEHICLE to MAP coordinates -> Calculate TransformedObservations
	//			Associate the nearest landmark with each TransformedObservations 
	//			Update Weight -> Apply Multivariate-Gaussian Probability to calculate P_xy
	//			Combine probabilities (multiply)

	double P_xy, scaleP_xy = 1./(2.*M_PI*std_landmark[0]*std_landmark[1]); //1.76839

	for (int i=0;i<num_particles;i++){ // for each PARTICLE

		std::vector<LandmarkObs> landmarks_in_sensor_range;
		std::vector<LandmarkObs> TransformedObservations;

		Particle this_particle = particles[i];

		double p_x = this_particle.x;
		double p_y = this_particle.y;
		double p_theta = this_particle.theta;

		// **
		// Get the list of landmarks within sensor_range of this particle: 
		// Using MAP COORDS
		// **
		for (int j=0;j<map_landmarks.landmark_list.size();j++){ // for each LANDMARK
			int lm_id = map_landmarks.landmark_list[j].id_i;
			double lm_x = map_landmarks.landmark_list[j].x_f;
			double lm_y = map_landmarks.landmark_list[j].y_f;
	
			if ( dist(p_x, p_y, lm_x, lm_y) < sensor_range) {
				LandmarkObs landmark = {lm_id, lm_x, lm_y};
				landmarks_in_sensor_range.push_back(landmark);
			}
		}

		for (int j=0;j<observations.size();j++){ // for each OBSERVATION

			LandmarkObs this_observation = observations[j];
			LandmarkObs tobs;

			// Transform the observation from VEHICLE to MAP COORDS
			tobs.x=this_observation.x*cos(p_theta) - this_observation.y*sin(p_theta) + p_x;
			tobs.y=this_observation.x*sin(p_theta) + this_observation.y*cos(p_theta) + p_y,
			tobs.id=this_observation.id;

			// And add it to the Transformed Observation vector
			TransformedObservations.push_back(tobs);

		} //for j OBSERVATION

//cout << "Particle: " << i << endl;
		// Having landmarks in range and observation coords in the same reference frame, make the nearest neighbor association
		dataAssociation(landmarks_in_sensor_range, TransformedObservations);

		// **
		// update the particle weight
		// **
		particles[i].weight = 1.;

		std::vector<int> associations;
		std::vector<double> sense_x;
		std::vector<double> sense_y;

		for (int j=0;j<TransformedObservations.size();j++){ // for each TRANSFORMED OBSERVATION

			LandmarkObs tobs = TransformedObservations[j];

			double mu_x=map_landmarks.landmark_list[tobs.id-1].x_f;
			double mu_y=map_landmarks.landmark_list[tobs.id-1].y_f;
//cout << tobs.x << "\t" << tobs.y << "\t" << mu_x << "\t" << mu_y << endl;
			//probability that this observation matches this landmark
			P_xy = scaleP_xy * exp(-(((tobs.x-mu_x)*(tobs.x-mu_x))/(2.*std_landmark[0]*std_landmark[0]) + ((tobs.y-mu_y)*(tobs.y-mu_y))/(2.*std_landmark[1]*std_landmark[1])));
			particles[i].weight *= P_xy; // Combine probabilities
//cout << "P_xy: " << P_xy << endl;
			associations.push_back(tobs.id);
			sense_x.push_back(tobs.x);
			sense_y.push_back(tobs.y);
		}

//cout << "particle weight:\t" << particles[i].weight << endl;

		particles[i] = SetAssociations(particles[i], associations, sense_x, sense_y);

		weights[i]=particles[i].weight;

	} /*for i (PARTICLE)*/

	for(int i=0;i<num_particles;i++){
		weights[i]=1/num_particles;
	}
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
