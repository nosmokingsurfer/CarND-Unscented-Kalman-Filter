#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"
#include "ground_truth_package.h"
#include "measurement_package.h"

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

  /************************************************************************/
  /* Scan the parameter space                                             */
  /************************************************************************/
  Eigen::VectorXd runSequence(const double nu_a_min, const double nu_a_max, const double nu_phidd_min, const double nu_phidd_max, const std::string input_file);


  /************************************************************************/
  /* Run one experiment                                                   */
  /************************************************************************/
  Eigen::VectorXd runOnce(const double nu_a, const double nu_phidd, std::vector<MeasurementPackage> &measurement_pack_list, std::vector<GroundTruthPackage> &gt_pack_list);

};

#endif /* TOOLS_H_ */
