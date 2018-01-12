#include <iostream>
#include "tools.h"
#include "ukf.h"


using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

  VectorXd result(4);
  result.fill(-1.0);

  if ((estimations.size() == 0) || (estimations.size() != ground_truth.size()))
    return result;

  for(int i = 0; i < static_cast<int>(estimations.size()); i++)
  {
    VectorXd estimation_temp = estimations[i];
    VectorXd ground_truth_temp = ground_truth[i];
    //std::cout << result << std::endl;
    //std::cout << estimation_temp << std::endl;
    //std::cout << ground_truth_temp << std::endl;
    result = result.array() + (estimation_temp - ground_truth_temp).array()*(estimation_temp - ground_truth_temp).array();
  }
  //std::cout << result << std::endl;

  result /= static_cast<int>(ground_truth.size());

  //std::cout << result << std::endl;

  result = result.array().sqrt();

  //std::cout << result << std::endl;
  return result;
}


VectorXd Tools::runSequence(const double nu_a_min, const double nu_a_max, const double nu_phidd_min, const double nu_phidd_max, const std::string input_file)
{
  VectorXd BestRMSE(4);
  BestRMSE.fill(0.0);

  double best_nu_a = -1;
  double best_nu_phidd = -1;
  double best_norm = 10;
  
  
  //Reading the data from file

  vector<MeasurementPackage> measurement_pack_list;
  vector<GroundTruthPackage> gt_pack_list;

  

  std::ifstream in_file_(input_file.c_str(), std::ifstream::in);
  std::string line;
  while (std::getline(in_file_, line))
  {
    std::string sensor_type;
    MeasurementPackage meas_package;
    GroundTruthPackage gt_package;
    std::istringstream iss(line);
    long long timestamp;

    // reads first element from the current line
    iss >> sensor_type;

    if (sensor_type.compare("L") == 0) {
      // laser measurement

      // read measurements at this timestamp
      meas_package.sensor_type_ = MeasurementPackage::LASER;
      meas_package.raw_measurements_ = VectorXd(2);
      float px;
      float py;
      iss >> px;
      iss >> py;
      meas_package.raw_measurements_ << px, py;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
    } else if (sensor_type.compare("R") == 0) {
      // radar measurement

      // read measurements at this timestamp
      meas_package.sensor_type_ = MeasurementPackage::RADAR;
      meas_package.raw_measurements_ = VectorXd(3);
      float ro;
      float phi;
      float ro_dot;
      iss >> ro;
      iss >> phi;
      iss >> ro_dot;
      meas_package.raw_measurements_ << ro, phi, ro_dot;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
    }

    // read ground truth data to compare later
    float x_gt;
    float y_gt;
    float vx_gt;
    float vy_gt;
    iss >> x_gt;
    iss >> y_gt;
    iss >> vx_gt;
    iss >> vy_gt;
    gt_package.gt_values_ = VectorXd(4);
    gt_package.gt_values_ << x_gt, y_gt, vx_gt, vy_gt;
    gt_pack_list.push_back(gt_package);
  }

  int N = 20;
  double nu_a_step = (nu_a_max - nu_a_min)/N;
  double nu_phidd_step = (nu_phidd_max - nu_phidd_min)/N;

  //Running the filter in cycle
  for(int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    {
      VectorXd tempResult;
      tempResult = runOnce(nu_a_min + i*nu_a_step, nu_phidd_min + j*nu_phidd_step, measurement_pack_list, gt_pack_list);
      double norm = tempResult.transpose()*tempResult;
      //std::cout << "norm = " << norm << std::endl;
      //if(norm < best_norm)
      //{
      //  best_norm = norm;
      //  best_nu_a = nu_a_min + i*nu_a_step;
      //  best_nu_phidd = nu_phidd_min + j*nu_phidd_step;
      //  BestRMSE = tempResult;
      //}

      if(tempResult[0] < 0.2)
      {
        if(tempResult[1] < 0.2)
        {
          if(tempResult[2] < 0.55)
          {
            if(tempResult[3] < 0.55)
            {
              best_norm = norm;
              best_nu_a = nu_a_min + i*nu_a_step;
              best_nu_phidd = nu_phidd_min + j*nu_phidd_step;
              BestRMSE = tempResult;      
            }
          }
        }
      }

    }
  }

  std::cout << "best nu_a = " << best_nu_a << std::endl;
  std::cout << "best nu phidd = " << best_nu_phidd << std::endl;
  std::cout << "BEST RMSE" << std::endl << BestRMSE << std::endl;
  std::cout << "best norm = " << best_norm << std::endl;

  return BestRMSE;
}

Eigen::VectorXd Tools::runOnce(const double nu_a, const double nu_phidd, vector<MeasurementPackage> &measurement_pack_list, vector<GroundTruthPackage> &gt_pack_list)
{
  VectorXd result(4);
  result.fill(0.0);

  UKF ukf;
  ukf.SetAccNoise(nu_a);
  ukf.SetPhiDDNoise(nu_phidd);
  vector<VectorXd> estimations;
  vector<VectorXd> ground_truth;


  size_t number_of_measurements = measurement_pack_list.size();

  for(size_t k = 0; k< number_of_measurements; k++)
  {
    
    ukf.ProcessMeasurement(measurement_pack_list[k]);
    VectorXd estimation(4);
    double x = ukf.x_[0];
    double y = ukf.x_[1];
    double vx = ukf.x_[2]*cos(ukf.x_[3]);
    double vy = ukf.x_[2]*sin(ukf.x_[3]);
    estimation << x, y, vx, vy;

    estimations.push_back(estimation);
    ground_truth.push_back(gt_pack_list[k].gt_values_);
  }
  
  

  Tools tools;
  result = tools.CalculateRMSE(estimations, ground_truth);
  std::cout << "nu_a = " << nu_a << std::endl;
  std::cout << "nu_phidd = " << nu_phidd << std::endl;
  std::cout << "result = " << std::endl << result << std::endl;
  
  return result;
}
