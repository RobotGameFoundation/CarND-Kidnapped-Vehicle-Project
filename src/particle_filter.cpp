/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  std::default_random_engine gen;
  num_particles = 100;  // TODO: Set the number of particles
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  for (int i = 0; i < num_particles; i++) {
    double sample_x = dist_x(gen);
    double sample_y = dist_y(gen);
    double sampel_theta = dist_theta(gen);
    Particle particle;
    particle.id = i;
    particle.x = sample_x;
    particle.y = sample_y;
    particle.theta = sampel_theta;
    particle.weight = 1;
    weights.push_back(particle.weight);
    particles.push_back(particle);
  }
  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  std::default_random_engine gen;
  // calculate particle x position will be inf if not limit yaw_rate  num
  if (fabs(yaw_rate)<0.0001) {
      yaw_rate = 0.0001;
    }
  
  for (auto& particle : particles) {
    double o_theta = particle.theta;
    particle.x += velocity/yaw_rate*(sin(o_theta+ yaw_rate*delta_t)-sin(o_theta));
    particle.y += velocity/yaw_rate*(cos(o_theta)-cos(o_theta+yaw_rate*delta_t));
    particle.theta += yaw_rate*delta_t;
    normal_distribution<double> dist_x(particle.x, std_pos[0]);
    normal_distribution<double> dist_y(particle.y, std_pos[1]);
    normal_distribution<double> dist_theta(particle.theta, std_pos[2]);
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    //std::cout << particle.x << std::endl;
  }
    

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, 
                                     std::vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

  
  for (auto& obj : observations) {
      double min_d = HUGE_VAL;
      for(auto& predict: predicted) { 
        //double distance = sqrt(pow(predicted.x-obj.x,2) + pow(predicted.y -obj.y,2));
        double distance = dist(obj.x, obj.y, predict.x, predict.y);
        if (distance < min_d) {
          min_d = distance;
          obj.id =  predict.id;
        }
      }
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
    
    
    double sig_x = std_landmark[0];
    double sig_y = std_landmark[1];
    double weight_sum = 0;
    for (auto& particle: particles) {
      
      double x = particle.x;
      double y = particle.y;
      double theta = particle.theta;

      vector<LandmarkObs> predicted;
      for (auto& landmark: map_landmarks.landmark_list) {
        
        int id = landmark.id_i;
      
        double x_landmark = landmark.x_f;
        double y_landmark = landmark.y_f;
        if (dist(x,y,x_landmark,y_landmark) < sensor_range) {
          LandmarkObs obs;
          obs.x = x_landmark;
          obs.y = y_landmark;
          obs.id = id;
          predicted.push_back(obs);
        }
       
      }
      vector<LandmarkObs> trans_obs;
      for (auto& obs: observations) {

        double x_obs = obs.x;
        double y_obs = obs.y;
        // transform to map x coordinate
        double x_map;
        x_map = x + (cos(theta) * x_obs) - (sin(theta) * y_obs);

        // transform to map y coordinate
        double y_map;
        y_map = y + (sin(theta) * x_obs) + (cos(theta) * y_obs);
        LandmarkObs p;
        p.x = x_map;
        p.y = y_map;
        trans_obs.push_back(p);
        
      }
      dataAssociation(predicted,trans_obs);
      particle.weight = 1.0;
      for (auto& obs: trans_obs) {
        
        double x_obs = obs.x;
        double y_obs = obs.y;
        int id_obs = obs.id;
        for (auto& pred: predicted) {
          double x_pred = pred.x;
          double y_pred = pred.y;
          int id_pred = pred.id;
          if (id_pred == id_obs) {
            particle.weight *= multiv_prob(sig_x, sig_y, x_obs, y_obs, x_pred, y_pred);
            
          }
        }
      }
      weight_sum += particle.weight;

    }

    
  for (int j=0; j<num_particles; j++) {
      Particle particle = particles[j];
      double weight = particle.weight/weight_sum;
      particle.weight = weight;
      weights[j] = weight;

    }

   
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
 
    
  std::vector<Particle> p3;
  double beta = 0.0;
  int index = rand()% num_particles;
  double max_weight = *max_element(weights.begin(), weights.end());
  for (int i=0; i<num_particles; i++) {
    double random = ((double) rand() / (RAND_MAX));
    beta += random*2.0*max_weight;
    while(beta > weights[index]) {
      beta -= weights[index];
      index = (index+1) % num_particles;
    }
    p3.push_back(particles[index]);
  }
  particles = p3;


}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}