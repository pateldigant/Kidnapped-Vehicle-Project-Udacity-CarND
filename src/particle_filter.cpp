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
#include <limits>
#include "map.h"
#include "particle_filter.h"
# define M_PI 3.14159265358979323846  /* pi */

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 100;

	default_random_engine gen;

	double std_x, std_y, std_theta; // Standard deviations for x, y, and theta
	std_x = std[0]; std_y = std[1]; std_theta = std[2];

	// This line creates a normal (Gaussian) distribution for x,y and theta.
	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);
	for(int i=0;i<num_particles;i++) {
		Particle particle;
		particle.id = i;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		particle.weight = 1.0;
		particles.push_back(particle);
		}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;

	double std_x, std_y, std_theta; // Standard deviations for x, y, and theta
	std_x = std_pos[0]; std_y = std_pos[1]; std_theta = std_pos[2];

	
	

	for(int i=0;i<num_particles;i++){
		if (fabs(yaw_rate) < 0.00001) {  
			particles[i].x += velocity * delta_t * cos(particles[i].theta);
			particles[i].y += velocity * delta_t * sin(particles[i].theta);
		}else{
			particles[i].x += velocity*( sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta)) / yaw_rate;
			particles[i].y += velocity*( cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t)) / yaw_rate; 
		}		
		particles[i].theta += yaw_rate*delta_t;
		// This line creates a normal (Gaussian) distribution for x,y and theta.
		normal_distribution<double> dist_x(particles[i].x, std_x);
		normal_distribution<double> dist_y(particles[i].y, std_y);
		normal_distribution<double> dist_theta(particles[i].theta, std_theta);

		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);

	}

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
	//std::vector<single_landmark_s> landmark_list = map_landmarks.landmark_list;
	
	for(int i=0;i<num_particles;i++){
			Particle p = particles[i];
			vector<int> landmarks;
			vector<double> sense_x,sense_y;
			double particle_weight = 1.0;

			for(int j=0;j<observations.size();j++){
				//convert the observations in map coordinate system
				double x_m =  p.x + (cos(p.theta)*observations[j].x) - (sin(p.theta)*observations[j].y);
				double y_m =  p.y + (sin(p.theta)*observations[j].x) + (cos(p.theta)*observations[j].y);
			
				sense_x.push_back(x_m);
				sense_y.push_back(y_m);

				//iterate over the landmarks to find the nearest one
				double best_distance = std::numeric_limits<double>::max();
				int best_index;
				double mu_x,mu_y,l_x,l_y,distance;
				for(int k=0;k<map_landmarks.landmark_list.size();k++){
					l_x = map_landmarks.landmark_list[k].x_f;
					 l_y = map_landmarks.landmark_list[k].y_f;
					 distance = (x_m - l_x)*(x_m - l_x) + (y_m - l_y)*(y_m - l_y);
					distance = sqrt(distance);
					if(distance<best_distance) best_distance=distance,best_index = map_landmarks.landmark_list[k].id_i,mu_x=l_x,mu_y=l_y;
				}
				landmarks.push_back(best_index);

				double gauss_norm = (1/(2 * M_PI * std_landmark[0] * std_landmark[1]));
				double t1,t2;
				
				t1 = x_m - mu_x;
				t2 = y_m - mu_y;
				double exponent= (pow(t1,2))/(2 * pow(std_landmark[0],2)) + (pow(t2,2))/(2 * std_landmark[1]*std_landmark[1]);

				double e_power = std::exp(-exponent);
				
				
				double weight= gauss_norm * e_power;
				
				particle_weight *= weight;
			}	
				p.weight = particle_weight;
				p.associations = landmarks;
				p.sense_x = sense_x;
				p.sense_y = sense_y;
				particles[i] = p;

		}

		



}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	double totalWeight = 0.0;
	for(int i=0;i<particles.size();i++) totalWeight += particles[i].weight;
	vector<double> alphas;
	double max_alpha = std::numeric_limits<double>::min();
	for(int i=0;i<particles.size();i++) alphas.push_back(particles[i].weight/totalWeight),max_alpha = max(max_alpha,particles[i].weight/totalWeight);
	
	std::default_random_engine generator;
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  	std::uniform_real_distribution<> distribution(0.0,2*max_alpha);
	std::uniform_int_distribution<> distr(0, num_particles-1);
	int index = distr(generator);
	double beta = 0.0;
	vector<Particle> tempParticles;
	for(int i=0;i<num_particles;i++){
		beta += distribution(gen); 
		while(alphas[index]< beta){
			beta -= alphas[index];
			index++;
			index %= num_particles;
		}
		tempParticles.push_back(particles[index]);
	}
	particles = tempParticles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
