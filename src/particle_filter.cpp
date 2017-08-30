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
#include <stdlib.h>

#include "particle_filter.h"

#define ZERO (0.000001F)
#define PI (3.141592F)

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
    num_particles = 50;

    particles.resize(num_particles);
    weights.resize(num_particles);

    default_random_engine gen;

    // Creates normal (Gaussian) distributions for x, y & theta.
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);

    for (int i = 0; i < num_particles; ++i) {
        // Sample from these normal distributions
        particles[i].x = dist_x(gen);
        particles[i].y = dist_y(gen);
        particles[i].theta = dist_theta(gen);
        particles[i].weight = 1.0f;
        weights[i] = 1.0f;
    }


#if 0  // DEBUG
    particles[0].x = 4;
    particles[0].y = 5;
    particles[0].theta = -PI/2;
    particles[0].weight = 1.0f;
    weights[0] = 1.0f;
#endif

    is_initialized = true;
}

void ParticleFilter::printParticles() {
    cout << "particles" << endl;
    cout << "---------" << endl;
    for (Particle p : particles) {
        cout << "  x:" << p.x << " y:" << p.y << " theta:" << p.theta << endl;
    }
    cout << endl;
}


void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // Update each particle's predicted pose (x, y, theta) based on the measured
    // controls (delta_t, velocity, yaw_rate) and add random Gaussian noise.
    // Using std::normal_distribution and std::default_random_engine for guassian noise:
    //   http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //   http://www.cplusplus.com/reference/random/default_random_engine/

    // Create normal (Gaussian) distributions for x, y & theta noise.
    normal_distribution<double> dist_x(0, std_pos[0]);
    normal_distribution<double> dist_y(0, std_pos[1]);
    normal_distribution<double> dist_theta(0, std_pos[2]);
    default_random_engine gen;  // random generator for the gaussian noise distribution.

    if (fabs(yaw_rate) > ZERO) {
        // for curved motion
        double delta_theta = yaw_rate * delta_t;    // pre-calculated for loop performance
        double huh = velocity / yaw_rate;           // pre-calculated for loop performance
        for (int i = 0; i < num_particles; ++i) {

            particles[i].x += huh * (sin(particles[i].theta + delta_theta) - sin(particles[i].theta));
            particles[i].x += dist_x(gen);          // add noise

            particles[i].y += huh * (cos(particles[i].theta) - cos(particles[i].theta + delta_theta));
            particles[i].y += dist_y(gen);          // add noise

            particles[i].theta += delta_theta;
            particles[i].theta += dist_theta(gen);  // add noise
        }
    } else {
        // for straight line motion (avoid divide by zero errors when yaw_rate = 0.
        double huh = velocity * delta_t;            // pre-calculated for loop performance
        for (int i = 0; i < num_particles; ++i) {
            particles[i].x += huh * cos(particles[i].theta);
            particles[i].x += dist_x(gen);          // add noise

            particles[i].y += huh * sin(particles[i].theta);
            particles[i].y += dist_y(gen);          // add noise
        }
    }
}


void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// For each observed landmark find the predicted landmark that is closest and associate the
	// observed landmark to it.
    //
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
    for (int j = 0; j < observations.size(); j++) {
        double nearest_dist = 10000.0f;
        int nearest = predicted.size() + 1;
        for (int i = 0; i < predicted.size(); i++) {
            double d = dist(predicted[i].x, predicted[i].y, observations[j].x, observations[j].y);
            if (d < nearest_dist) {
                nearest_dist = d;
                nearest = i;
            }
        }
        if (nearest < predicted.size()) {
            observations[j].id = predicted[nearest].id;  // this is the index into the landmark list in the map.
        }
    }
}

// Build a list of landmarks that are within sensor range of the particle.
void ParticleFilter::predict_landmarks(Particle particle, double sensor_range, Map map, std::vector<LandmarkObs>& predicted_landmarks) {

    predicted_landmarks.clear();

    for (int i = 0; i < map.landmark_list.size(); i++) {
        // distance from particle to map landmark
        double distance = dist(particle.x, particle.y, map.landmark_list[i].x_f, map.landmark_list[i].y_f);
        // include all landmarks that are within sensor range of the particle
//        if (distance < sensor_range) {
            // add landmark to predicted landmarks
            LandmarkObs pred;
            pred.x = map.landmark_list[i].x_f;
            pred.y = map.landmark_list[i].y_f;
            pred.id = i;    // save index to landmark list in the map.
            predicted_landmarks.push_back(pred);
//        }
    }
}

// Transform an observation from car's perspective (relative to car) to particle's perspective (relative to map).
LandmarkObs transform_coords(Particle particle, LandmarkObs observation) {
#if 0 // DEBUG
    cout << "  observation x: " << observation.x << " y: " << observation.y << endl;
#endif
    LandmarkObs transform;
    transform.x = particle.x + (observation.x * cos(particle.theta)) - (observation.y * sin(particle.theta));
//    if (transform.x < ZERO)
//        transform.x = 0.0f;
    transform.y = particle.y + observation.x * sin(particle.theta) + observation.y * cos(particle.theta);
//    if (transform.y < ZERO)
//        transform.y = 0.0f;
#if 0 // DEBUG
    cout << "  transform   x: " << transform.x << " y: " << transform.y << endl << endl;
#endif

    return(transform);
}


double measurement_prob(LandmarkObs observation, Map::single_landmark_s landmark, double std_landmark[]) {
    double err_x = observation.x - landmark.x_f;
    double err_y = observation.y - landmark.y_f;
    double p1 = (err_x * err_x) / (2.0f * std_landmark[0] * std_landmark[0]) + (err_y * err_y) / (2.0f * std_landmark[1] * std_landmark[1]);
    double prob = exp(-p1) / (2.0f * PI * std_landmark[0] * std_landmark[1]);
    return (prob);
}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   std::vector<LandmarkObs> observations, Map map_landmarks) {
    // Update the weights of each particle using a multi-variate Gaussian distribution. You can read
	  // more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    //
    // The observations are given in the VEHICLE'S coordinate system. The particles are located
    // according to the MAP'S coordinate system. It's necessary to transform between the two systems.
    // This transformation must support both rotation and translation.
    //
    // The following is a good resource for the theory:
    //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
    // and the following is a good resource for the actual equation to implement (look at equation 3.33):
    //   http://planning.cs.uiuc.edu/node99.html
    for (int i = 0; i < num_particles; i++) {
        // Transform the observations from the car's perspective to the particle's perspective.
        // The transform calculates where the observation would appear on the map if observed
        // from the particle's perspective.

#if 0
        cout << "transformed observations for particle " << i << endl;
        cout << "-------------------------------------" << endl;
        cout << "  position    x: " << particles[i].x << " y: " << particles[i].y << endl << endl;
#endif
        std::vector<LandmarkObs> transformed_observations(observations.size());
        for (int j = 0 ; j < observations.size(); j++) {
            transformed_observations[j] = transform_coords(particles[i], observations[j]);
        }

        // Generate list of all landmarks that are within sensor range of the particle's predicted position
        std::vector<LandmarkObs> predicted_landmarks;
        predict_landmarks(particles[i], sensor_range, map_landmarks, predicted_landmarks);

        // Find associations between the transformed observations and the particle's predicted landmarks.
        dataAssociation(predicted_landmarks, transformed_observations);

        // Calculate error between the particle's transformed observations and the associated map landmarks.
        particles[i].weight = 1.0f;
#if 0
        cout << "particle(" << i << ")" << endl;
        cout << "------------------" << endl;
        int k = 0;
#endif
        for (LandmarkObs observation : transformed_observations) {
            // calculate the multivariate gaussian probability
            double p = measurement_prob(observation, map_landmarks.landmark_list[observation.id], std_landmark);
            particles[i].weight *= p;
#if 0
            cout << "  observation(" << k++ << ") probability: " << p << endl;
#endif
        }
        weights[i] = particles[i].weight;
#if 0
        cout << endl << "  weight: " << weights[i] << endl << endl;
#endif
    }
}

void ParticleFilter::resample() {
    // Resample particles with replacement with probability proportional to their weight.
#if 0
    // resample using resampling wheel
    float fraction = ((float)rand()/(float)(RAND_MAX));
    int index = int(fraction * (float)num_particles);  // pick a random starting point position on the wheel
    double beta = 0.0;
    // resampling wheel is weighted toward selecting the highest measurement probability
    double mw = *max_element(weights.begin(), weights.end());
    cout << "max element: " << mw << endl;
    cout << "resampling" << endl;
    cout << "----------" << endl;
    std::vector<Particle> resampled;
    for (int i = 0; i < num_particles; i++) {
        fraction = ((float)rand()/(float)(RAND_MAX)); // generates a float between 0.0 and 1.0
        beta += fraction * 2.0 * mw;
        while (beta > weights[index]) {
            beta -= weights[index];
            index = (index + 1) % num_particles;
        }
        cout << "  selected " << index << endl;
        resampled.push_back(particles[index]);
    }
#else
    // resample using std::discrete_distribution
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> distribution(weights.begin(), weights.end());
    std::vector<Particle> resampled;
    int index;
    for(int n=0; n<num_particles; n++) {
        index = distribution(gen);
        resampled.push_back(particles[index]);
    }
#endif

    particles.clear();
    for (Particle p : resampled) {
      particles.push_back(p);
    }
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
